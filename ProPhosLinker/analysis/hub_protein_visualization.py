#!/usr/bin/env python3
"""
ProPhosLinker: Hub protein expansion and visualization module for the ProPhosLinker pipeline.

This module provides the `hub_protein_expansion()` function (part of ProPhosLinker toolkit), 
which identifies top hub proteins from knowledge graph functional networks and generates
individualized expanded network visualizations for each hub.

Core responsibilities:
1) **Multi-metric hub identification**: Computes composite centrality scores using:
   * Degree centrality (40% weight)
   * Betweenness centrality (30% weight) 
   * Closeness centrality (20% weight)
   * Eigenvector centrality (10% weight, with robust fallback handling)

2) **Robust eigenvector centrality**: Handles convergence failures and disconnected graphs
   by falling back to component-wise computation or degree centrality.

3) **Per-hub network expansion**: For each top-N hub protein:
   * Creates dedicated output directory (`Hub_Protein_<HUBNAME>`)
   * Calls R visualization script with comprehensive styling parameters
   * Supports network modes (ALL), layouts (kk), and phospho-site limits

4) **Rich visualization parameters**:
   * Node colors by regulation state (up/down/nonsig/notdet)
   * Disease-contextualized border colors
   * 14 edge-type specific colors (association, activation, catalysis, etc.)
   * Function enrichment color gradients

Typical usage:
- Called from `FunctionalAnalysis.hub_protein_functional_network()` after knowledge
  graph network construction (step 3.3).
- Expects `Node_diff_info.tsv` and `Edge_diff_info.tsv` from previous step.
"""


from pathlib import Path
from typing import Dict
import pandas as pd
import networkx as nx
import subprocess 
import os
import sys

def hub_protein_expansion(node_path,edge_path,
    outdir_parent,
    top_n=6,
    d_num = 5,
    network_mode = "ALL",
    network_layout = 'kk',
    omics1_name = 'Pro',
    omics2_name ='Phos',
    identified_type = 'UNIPROT',
    r_script_path = "E:/pro_phos_V2/scripts/3.2network_hubgen_visualization.R",
    function_enrich_rscript_path = '',
    hub_protein_network_node_up = "#B83A2D",
    hub_protein_network_node_down = "#657f68",
    hub_protein_network_node_nonsig = "#E3C79F",
    hub_protein_network_node_notdet = "grey90",
    hub_protein_network_node_disease_border_color = "#535c54",
    hub_protein_network_node_notdisease_border_color = "#acabab",
    hub_protein_network_function_enrich_color_gradient_low = "#175663",
    hub_protein_network_function_enrich_color_gradient_high = "#90362d",
    # Edge type color mapping (used for visualization)
    hub_protein_network_edge_type_colors = {
        "association": "#d7d6d6",           # light gray
        "physical association": "#838181",  # light gray (same as binding)
        "binding": "#838181",               # light gray
        "direct interaction": "#d7d6d6",    # light gray
        
        "activation": "#d62c0c",            # orange-red
        "catalysis": "#7f00ff",             # light purple
        "proximity": "#FF9301",             # orange
        "reaction": "#114335",              # dark green
        "phosphorylation": "#2FBE95",       # teal
        "dephosphorylation": "#6B8E23",     # olive
        
        "ptmod": "#8C97D6",                 # light purple
        "inhibition": "#0cb6d6",            # deep blue
        "expression": "#FCF402",            # yellow
        "regulation": "#e4dadd",            # mauve gray
        
        "colocalization": "#4c95cd",        # sky blue
        "covalent binding": "#716F74",      # violet
        "ubiquitination": "#FF4500",        # orange-red
        "multiRel": "#9a6728"               # brown
    }
    ):
    try:
        # Create the output directory if it does not exist
        if not outdir_parent.exists():
            outdir_parent.mkdir(parents=True, exist_ok=True)
        try:
            outdir_parent.mkdir(parents=True, exist_ok=True)
            # Status message: output directory is ready
            print(f"✅ Output directory ready: {outdir_parent}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {outdir_parent}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")
        node_df = pd.read_csv(node_path,sep='\t')
        edge_df = pd.read_csv(edge_path,sep='\t')
       
        G = nx.Graph()
        for _, row in node_df.iterrows():
            G.add_node(row['NodeName'])
        for _, row in edge_df.iterrows():
            G.add_edge(row['Source'], row['Target'])


        # Compute multiple centrality metrics
        def get_top_hub_nodes(G, top_n=6):
            """Select hub nodes using a composite score across multiple centrality measures.

            The implementation is robust to eigenvector centrality convergence issues.
            """
            # 1. Degree centrality
            degree_centrality = nx.degree_centrality(G)
            # 2. Betweenness centrality
            betweenness_centrality = nx.betweenness_centrality(G)
            # 3. Closeness centrality
            closeness_centrality = nx.closeness_centrality(G)

            # 4. Eigenvector centrality (robust implementation)
            eigenvector_centrality = {n: 0.0 for n in G.nodes()}
            try:
                # First try power iteration with increased max iterations and tighter tolerance
                eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000, tol=1.0e-06)
            except nx.PowerIterationFailedConvergence:
                try:
                    # If power iteration fails, try numpy-based eigendecomposition
                    if nx.is_connected(G):
                        eigenvector_centrality = nx.eigenvector_centrality_numpy(G)
                    else:
                        # For disconnected graphs, compute per component and merge
                        eigenvector_centrality = {}
                        for comp in nx.connected_components(G):
                            sub = G.subgraph(comp)
                            if len(sub) == 1:
                                # For isolated node, assign 0.0 (or 1.0 if desired)
                                node = next(iter(sub.nodes()))
                                eigenvector_centrality[node] = 0.0
                                continue
                            try:
                                sub_ev = nx.eigenvector_centrality_numpy(sub)
                                eigenvector_centrality.update(sub_ev)
                            except Exception:
                                # Fallback to degree centrality for this component to avoid breaking the workflow
                                for n in sub.nodes():
                                    eigenvector_centrality[n] = degree_centrality.get(n, 0.0)
                except Exception:
                    # If all methods fail, fall back to degree centrality as a last resort
                    print("Warning: eigenvector centrality calculation failed. Falling back to degree centrality.", file=sys.stderr)
                    eigenvector_centrality = {n: degree_centrality.get(n, 0.0) for n in G.nodes()}

            # Composite score: weighted average (preserving original weights)
            nodes = list(G.nodes())
            composite_scores = {}
            for node in nodes:
                composite_score = (
                    degree_centrality.get(node, 0.0) * 0.4 +
                    betweenness_centrality.get(node, 0.0) * 0.3 +
                    closeness_centrality.get(node, 0.0) * 0.2 +
                    eigenvector_centrality.get(node, 0.0) * 0.1
                )
                composite_scores[node] = composite_score

            # Select top `top_n` hub nodes by composite score
            top_hub_nodes = sorted(composite_scores.items(), key=lambda x: x[1], reverse=True)[:top_n]
            return top_hub_nodes, {
                'degree': degree_centrality,
                'betweenness': betweenness_centrality,
                'closeness': closeness_centrality,
                'eigenvector': eigenvector_centrality,
                'composite': composite_scores
            }

        top6_hub_nodes, centrality_metrics = get_top_hub_nodes(G, top_n)

        for hub_protein in top6_hub_nodes:
            outdir = outdir_parent / Path("Hub_Protein_"+hub_protein[0])
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            d_num = 5
            network_mode = "ALL"
            network_layout = 'kk'
            # network_layout = 'tree'
            visualization_SEEDNODEID = hub_protein[0]
            def to_r_path(path):
                return str(path).replace("\\", "/")
            def edge_colors_to_param(d: Dict[str, str], entry_sep: str = ';', kv_sep: str = ':') -> str:
                parts = [f"{k.strip()}{kv_sep}{v.strip()}" for k, v in d.items()]
                return entry_sep.join(parts)
            
            cmd = [
                "Rscript",
                to_r_path(r_script_path),
                "--nodes", to_r_path(node_path),
                "--edges", to_r_path(edge_path),
                "--max_phosphoSite_displayed", str(d_num),
                "--outdir", to_r_path(outdir),
                "--network_mode", network_mode,
                "--SEEDNODEID", visualization_SEEDNODEID,
                "--network_layout", network_layout,
                "--omics1_name", omics1_name,
                "--omics2_name", omics2_name,
                "--enrich_fromType", identified_type,
                "--node_up", hub_protein_network_node_up,
                "--node_down", hub_protein_network_node_down,
                "--node_nonsig", hub_protein_network_node_nonsig,
                "--node_notdet", hub_protein_network_node_notdet,
                "--node_disease_border_color", hub_protein_network_node_disease_border_color,
                "--node_notdisease_border_color", hub_protein_network_node_notdisease_border_color,
                "--function_enrich_color_gradient_low", hub_protein_network_function_enrich_color_gradient_low,
                "--function_enrich_color_gradient_high", hub_protein_network_function_enrich_color_gradient_high,
                "--edge_type_colors", edge_colors_to_param(hub_protein_network_edge_type_colors),
                "--function_enrich_rscript_path", to_r_path(function_enrich_rscript_path)
            ]
           
            full_cmd = cmd
            
            print(f"📋 Running command:\n{' '.join(full_cmd)}")
            result = subprocess.run(full_cmd, capture_output=True, text=True)

            print("STDOUT:\n", result.stdout)
            print("STDERR:\n", result.stderr)
        return True
    except Exception as e:
        print(f"hub_protein_expansion error: {e}", file=sys.stderr)
        return False
