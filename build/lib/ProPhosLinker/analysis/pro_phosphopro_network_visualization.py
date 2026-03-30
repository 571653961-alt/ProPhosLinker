#!/usr/bin/env python3
"""
ProPhosLinker: Network Visualization Module (3.4)

This module generates publication-ready visualizations of protein-phosphoprotein 
functional interaction networks from knowledge graph analysis results. It provides:

Key Features:
- Module-specific subnetwork extraction and filtering
- Multi-mode network filtering (KS, LR, TF, DEP, SEEDNODE)
- Centrality-based node ranking (9 algorithms: betweenness, degree, closeness, etc.)
- Phosphosite aggregation with configurable display limits
- Integration with R visualization engine via subprocess execution
- Automated parameter passing for consistent styling
- Support for disease association and differential expression visualization

Workflow:
1. Module extraction → Node filtering → Centrality ranking → Subgraph creation
2. TSV export for R visualization → Dynamic R script execution
3. Comprehensive styling with edge type colors and node state encoding

Author: ProPhosLinker Team
License: MIT
Dependencies: networkx, pandas, subprocess
"""

import os
import sys
from pathlib import Path
from typing import Union, Dict, Literal
import pandas as pd
import re
import networkx as nx
import subprocess


class NetworkConfig:
    """
    Comprehensive configuration class for network visualization pipeline.
    
    Centralizes all parameters for network extraction, filtering, centrality 
    analysis, and R visualization integration.
    """
    
    def __init__(
        self,
        *,
        # Input files
        profile: Union[str, Path, pd.DataFrame],
        phosphoprofile: Union[str, Path, pd.DataFrame],
        prodiff: Union[str, Path, pd.DataFrame],
        phosphoprodiff: Union[str, Path, pd.DataFrame],
        phosphopro_pro_cor: Union[str, Path, Dict],
        # Output directory
        outdir: Union[str, Path],
        # Network data parameters
        nodefile: Union[str, Path, pd.DataFrame] = None,
        edgefile: Union[str, Path, pd.DataFrame] = None,
        module_info: Union[str, Path, pd.DataFrame] = None,
        # Omics names
        omics1_name: str = 'Pro',
        omics2_name: str = 'Phos',
        # Network construction parameters
        FC: float = 1.5,
        disease: str = '',
        kgbased_network_top_n: int = 6,
        # Network analysis parameters
        analysis_network_mode: Literal["ALL", "KS", "LR", "TF", 
                                     "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'DEP',
        analysis_SEEDNODEID: str = None,
        comm_detection: Literal["louvain", "leiden", "infomap", 
                               "fastgreedy", "walktrap", "LPA"] = 'infomap',
        # Visualization parameters
        top_nodes_visualization_num: int = 20,
        module_id: int = 0,
        max_phosphoSite_displayed: int = 5,
        visualization_network_mode: Literal["ALL", "KS", "LR", "TF", 
                                          "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'ALL',
        visualization_SEEDNODEID: str = None,
        node_filtering: Literal["betweenness", "degree", "closeness", "harmonic",
                               "eigenvector", "pagerank", "alpha", "hub", "authority"] = 'betweenness',
        network_layout: Literal["fr", "kk", "dh", "stress", "tree", "gem", 
                               "graphopt", "lgl", "circle", "grid"] = 'kk',
        oudir_default_name: str = "3.2KGbased_Functional_Network",
        # Node colors (UP/DOWN/Non-sig/Not-detected)
        node_up: str = "#B83A2D",                    # Upregulated red
        node_down: str = "#657f68",                  # Downregulated green-gray
        node_nonsig: str = "#E3C79F",                # Non-significant yellow
        node_notdet: str = "grey90",                 # Not detected light gray
        node_disease_border_color: str = "#535c54",  # Disease-related border
        node_notdisease_border_color: str = "#acabab", # Non-disease border
        function_enrich_color_gradient_low: str = "#175663",   # Low enrichment blue
        function_enrich_color_gradient_high: str = "#90362d",  # High enrichment red
        # Edge type colors
        edge_type_colors: Dict[str, str] = {
            "association": "#d7d6d6",           # Light gray
            "physical association": "#838181",  # Medium gray
            "binding": "#838181",               # Medium gray
            "direct interaction": "#d7d6d6",    # Light gray
            
            "activation": "#d62c0c",            # Red-orange
            "catalysis": "#7f00ff",             # Purple
            "proximity": "#FF9301",             # Orange
            "reaction": "#114335",              # Dark green
            "phosphorylation": "#2FBE95",       # Teal
            "dephosphorylation": "#6B8E23",     # Olive green
            
            "ptmod": "#8C97D6",                 # Light purple
            "inhibition": "#0cb6d6",            # Cyan
            "expression": "#FCF402",            # Yellow
            "regulation": "#e4dadd",            # Purple-gray
            
            "colocalization": "#4c95cd",        # Sky blue
            "covalent binding": "#716F74",      # Violet
            "ubiquitination": "#FF4500",        # Orange-red
            "multiRel": "#9a6728"               # Brown
        },
        # R script paths
        script_path: str = './scripts',
        kgbased_functional_network_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_visualization.R',
        kgbased_functional_network_community_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_community_visualisation.R',
        function_enrich_rscript_path: Union[str, Path] = '3.2functional_enrichment_function.R',
        # Neo4j connection
        uri: str = "bolt://localhost:7687",
        username: str = "neo4j",
        password: str = "neo4j"
    ):
        # Store all configuration parameters
        for key, value in locals().items():
            if key != 'self':
                setattr(self, key, value)
        
        # Initialize and validate
        self._validate_params()
        self._data_load()
    
    def _validate_params(self):
        """Validate inputs and setup output directory structure."""
        # Convert outdir to Path and validate
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        elif self.outdir is None:
            self.outdir = Path(os.getcwd())
        elif not isinstance(self.outdir, Path):
            raise TypeError("outdir must be string or Path object")
        
        # Validate output directory permissions
        if self.outdir.exists():
            if self.outdir.is_file():
                raise ValueError(f"Output path is file, not directory: {self.outdir}")
            if not os.access(self.outdir, os.W_OK):
                raise ValueError(f"Output directory not writable: {self.outdir}")
        
        # Create default output subdirectory
        self.outdir = self.outdir / Path(self.oudir_default_name)
        output_dir = self.outdir
        
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f"✅ Output directory ready: {output_dir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {output_dir}")
        except Exception as e:
            raise RuntimeError(f"Directory creation failed: {str(e)}")
        
        # Validate input files
        for attr in ['profile', 'phosphoprofile', 'prodiff', 'phosphoprodiff', 'phosphopro_pro_cor']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Path(value))
            elif not isinstance(value, (Path, pd.DataFrame)):
                raise TypeError(f"{attr} must be string, Path, or DataFrame")
            
            # File validation
            if isinstance(value, Path):
                if not value.exists():
                    raise FileNotFoundError(f"{value} does not exist")
                if not value.is_file():
                    raise ValueError(f"{value} is not valid file")
            
            # DataFrame validation
            if isinstance(value, pd.DataFrame):
                if value.empty:
                    raise ValueError(f"{attr} DataFrame is empty")
                if not all(col in value.columns for col in ['ID', 'Name']):
                    raise ValueError(f"{attr} must contain 'ID' and 'Name' columns")
        
        # Validate fold-change threshold
        if self.FC <= 0:
            raise ValueError("Fold change (log2FC) must be > 0")
    
    def _data_load(self):
        """Load and preprocess all input data files."""
        def datafile_load_validation(file):
            """Load expression data with numeric conversion."""
            if not isinstance(file, pd.DataFrame):
                df = pd.read_csv(file, sep='\t', index_col=0, 
                               na_values=["NA", "-", "null", ""])
                df = df.apply(pd.to_numeric, errors='coerce')
            else:
                df = file
            
            if df.empty:
                raise ValueError(f"{file} expression data is empty")
            if df.index.duplicated().any():
                raise ValueError(f"{file} contains duplicate indices")
            return df
        
        def difffile_load_validation(file):
            """Load differential expression results."""
            if not isinstance(file, pd.DataFrame):
                df = pd.read_csv(file, sep='\t', index_col=0, 
                               na_values=["NA", "-", "null", ""])
            else:
                df = file
            return df
        
        # Load all datasets
        self.profile = datafile_load_validation(self.profile)
        self.phosphoprofile = datafile_load_validation(self.phosphoprofile)
        self.prodiff = difffile_load_validation(self.prodiff)
        self.phosphoprodiff = difffile_load_validation(self.phosphoprodiff)
        
        # Load phosphosite-to-protein mapping
        if not isinstance(self.phosphopro_pro_cor, dict):
            self.phosphopro_pro_cor = pd.read_csv(
                self.phosphopro_pro_cor, sep='\t', header=None,
                names=['value', 'key']
            ).set_index('key')['value'].to_dict()


def subnetwork_construction(config: NetworkConfig):
    """
    Construct visualization-ready subnetwork for specific module.
    
    Pipeline:
    1. Extract module nodes
    2. Filter phosphosites by mode and display limit
    3. Build directed graph from filtered edges
    4. Apply network mode filtering (KS/LR/TF/SEEDNODE)
    5. Rank nodes by centrality → Select top N
    6. Export filtered nodes/edges for R visualization
    
    Args:
        config: NetworkConfig with analysis parameters
        
    Returns:
        tuple: (nodes_path, edges_path) TSV file paths
    """
    # Setup output paths
    viz_nodes_path = config.outdir / 'nodes_visualization.tsv'
    viz_edges_path = config.outdir / 'edges_visualization.tsv'
    
    # Load input data
    node_raw = config.nodefile
    edge_raw = config.edgefile
    module_info = config.module_info
    module_id = config.module_id
    top_n_nodes = config.top_nodes_visualization_num
    network_mode = config.visualization_network_mode
    seed_node = config.visualization_SEEDNODEID
    max_phosphosites = config.max_phosphoSite_displayed
    centrality_metric = config.node_filtering
    
    # Extract module nodes
    module_items = module_info.loc[module_info['ModuleID'] == module_id, 'ModuleItems'].iloc[0]
    module_nodes = module_items.split(';')
    
    # Validate SEEDNODE mode
    if network_mode == 'SEEDNODE':
        if not seed_node:
            raise ValueError("SEEDNODEID cannot be empty in SEEDNODE mode")
        if seed_node not in module_nodes:
            raise ValueError(f"SEEDNODEID '{seed_node}' not in module {module_id}")
    
    # Filter nodes to module
    module_node_data = node_raw[node_raw['NodeName'].isin(module_nodes)].copy()
    
    # 1. Phosphosite filtering by network mode
    def filter_nodes_by_mode(node_df, mode, max_phos):
        """Filter and rank phosphosites within display limits."""
        # Ensure numeric types
        node_df['Phosphopro_FC'] = pd.to_numeric(node_df['Phosphopro_FC'], errors='coerce')
        
        # Mode-specific filtering
        if mode == 'DEP':
            mask = (node_df['Pro_class'].isin(['UP', 'DOWN'])) & \
                   (node_df['Phosphopro_class'].isin(['UP', 'DOWN']))
        elif mode == 'DEPUP':
            mask = (node_df['Pro_class'] == 'UP') & (node_df['Phosphopro_class'] == 'UP')
        elif mode == 'DEPDOWN':
            mask = (node_df['Pro_class'] == 'DOWN') & (node_df['Phosphopro_class'] == 'DOWN')
        else:
            mask = pd.Series(True, index=node_df.index)
        
        filtered = node_df[mask].copy()
        
        # Rank phosphosites by absolute fold-change, limit per protein
        result = filtered.sort_values(['NodeName', 'Phosphopro_FC'], 
                                    key=lambda x: abs(x) if x.name == 'Phosphopro_FC' else x,
                                    ascending=[True, False]).groupby('NodeName', group_keys=False).head(max_phos)
        
        # Add ranking index
        result['phosphoSite_filter_index'] = result.groupby('NodeName').cumcount() + 1
        
        if result.empty:
            raise ValueError(f"No nodes match filtering criteria for mode '{mode}'")
        
        return result
    
    filtered_nodes = filter_nodes_by_mode(module_node_data, network_mode, max_phosphosites)
    
    # 2. Construct initial graph
    nodes = filtered_nodes['NodeName'].unique()
    module_edges = edge_raw[
        edge_raw['Source'].isin(nodes) & edge_raw['Target'].isin(nodes)
    ].copy()
    
    # Create directed graph
    G = nx.from_pandas_edgelist(
        module_edges, source='Source', target='Target',
        edge_attr=['Interaction_Type'], create_using=nx.DiGraph()
    )
    print(f"✅ Initial graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges", file=sys.stderr)
    
    # 3. Network mode edge filtering
    def filter_edges_by_mode(G, mode, seed_id=None):
        """Apply mode-specific edge filtering."""
        if mode == 'KS':
            # Kinase-substrate (phosphorylation)
            edges_keep = [(u, v) for u, v, d in G.edges(data=True)
                         if 'Interaction_Type' in d and re.search('phospho', d['Interaction_Type'], re.IGNORECASE)]
            G_filtered = G.edge_subgraph(edges_keep).copy()
        
        elif mode == 'LR':
            # Ligand-receptor (binding)
            edges_keep = [(u, v) for u, v, d in G.edges(data=True)
                         if 'Interaction_Type' in d and re.search('binding', d['Interaction_Type'], re.IGNORECASE)]
            G_filtered = G.edge_subgraph(edges_keep).copy()
        
        elif mode == 'TF':
            # Transcription factor interactions
            tf_patterns = 'activation|regulation|repression|inhibition'
            edges_keep = [(u, v) for u, v, d in G.edges(data=True)
                         if 'Interaction_Type' in d and re.search(tf_patterns, d['Interaction_Type'], re.IGNORECASE)]
            G_filtered = G.edge_subgraph(edges_keep).copy()
        
        elif mode == 'SEEDNODE' and seed_id:
            # Seed node neighborhood
            neighbors1 = set(G.neighbors(seed_id))
            neighbors1.add(seed_id)
            
            if len(neighbors1) < 10:
                neighbors2 = set()
                for node in neighbors1:
                    neighbors2.update(G.neighbors(node))
                sub_nodes = neighbors1 | neighbors2
            else:
                sub_nodes = neighbors1
            
            G_filtered = G.subgraph(sub_nodes).copy()
        else:
            G_filtered = G.copy()
        
        # Remove isolates
        G_filtered.remove_nodes_from(list(nx.isolates(G_filtered)))
        return G_filtered
    
    G_filtered = filter_edges_by_mode(G, network_mode, seed_node)
    print(f"✅ Filtered graph: {G_filtered.number_of_nodes()} nodes, {G_filtered.number_of_edges()} edges", file=sys.stderr)
    
    # 4. Centrality-based node ranking
    def compute_centrality_top_nodes(G, metric, top_n):
        """Calculate centrality and return top N nodes."""
        centrality_funcs = {
            'betweenness': lambda g: nx.betweenness_centrality(g, normalized=False),
            'degree': nx.degree_centrality,
            'closeness': nx.closeness_centrality,
            'harmonic': nx.harmonic_centrality,
            'eigenvector': lambda g: nx.eigenvector_centrality(g, max_iter=1000),
            'pagerank': nx.pagerank,
            'alpha': lambda g: nx.katz_centrality(g, alpha=0.85),
            'hub': lambda g: nx.hits(g)[0],
            'authority': lambda g: nx.hits(g)[1]
        }
        
        if metric not in centrality_funcs:
            raise ValueError(f"Unknown centrality metric: {metric}")
        
        centrality = centrality_funcs[metric](G)
        top_nodes = [node for node, _ in sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:top_n]]
        
        if not top_nodes:
            raise ValueError("No nodes found after centrality filtering")
        
        return top_nodes
    
    def create_centrality_subgraph(G, metric, top_n, mode):
        """Create subgraph with top centrality nodes."""
        try:
            top_node_ids = compute_centrality_top_nodes(G, metric, top_n)
            subG = G.subgraph(top_node_ids).copy()
            
            # Remove isolates except for differential expression modes
            if mode not in ["DEP", "DEPUP", "DEPDOWN"]:
                isolates = [n for n in subG.nodes if subG.degree(n) == 0]
                subG.remove_nodes_from(isolates)
            
            return subG
        except Exception as e:
            print(f"❌ Subgraph creation failed: {str(e)}", file=sys.stderr)
            return nx.DiGraph()
    
    # Create final visualization subgraph
    viz_subgraph = create_centrality_subgraph(
        G_filtered, centrality_metric, top_n_nodes, network_mode
    )
    print(f"✅ Visualization subgraph: {viz_subgraph.number_of_nodes()} nodes, {viz_subgraph.number_of_edges()} edges", file=sys.stderr)
    
    # 5. Export for R visualization
    nx.to_pandas_edgelist(viz_subgraph).to_csv(
        viz_edges_path, sep='\t', index=False, encoding='utf-8'
    )
    viz_nodes = config.nodefile[
        config.nodefile['NodeName'].isin(viz_subgraph.nodes)
    ].to_csv(viz_nodes_path, sep='\t', index=False, encoding='utf-8')
    
    return viz_nodes_path, viz_edges_path


def NetworkVisualization(config: NetworkConfig):
    """
    Execute complete network visualization pipeline for top modules.
    
    1. Creates module-specific visualization directories
    2. Constructs optimized subnetworks for each module  
    3. Launches R visualization scripts with full parameter set
    4. Handles edge type color mapping and path conversion
    
    Args:
        config: NetworkConfig with visualization parameters
        
    Side Effects:
        - Creates module-specific output directories
        - Generates TSV files for R processing
        - Executes R scripts for publication-ready figures
    """
    # Setup visualization output directory
    viz_dir = config.outdir / f"{config.comm_detection}_community_detection_modules_{config.module_id}_visualization"
    viz_dir.mkdir(parents=True, exist_ok=True)
    config.outdir = viz_dir
    
    network_mode = config.visualization_network_mode
    seed_node = config.visualization_SEEDNODEID or 'None'
    centrality_metric = config.node_filtering
    top_nodes = config.top_nodes_visualization_num
    max_phos_display = config.max_phosphoSite_displayed
    layout_algo = config.network_layout
    
    def to_r_path(path: Path) -> str:
        """Convert paths to R-compatible format."""
        return str(path).replace("\\", "/")
    
    def edge_colors_to_param(colors: Dict[str, str], entry_sep: str = ';', kv_sep: str = ':') -> str:
        """Format edge type colors for R parameter passing."""
        parts = [f"{k.strip()}{kv_sep}{v.strip()}" for k, v in colors.items()]
        return entry_sep.join(parts)
    
    # Process top N modules only
    module_ids = config.module_info['ModuleID'].tolist()
    if len(module_ids) > config.kgbased_network_top_n:
        module_ids = module_ids[:config.kgbased_network_top_n]
    
    for module_id in module_ids:
        print(f"🎨 Visualizing module {module_id}...", file=sys.stderr)
        
        # Update config for current module
        config.module_id = module_id
        
        # Construct module-specific subnetwork
        node_viz_path, edge_viz_path = subnetwork_construction(config)
        
        # Build comprehensive R command
        r_cmd = [
            "Rscript",
            to_r_path(config.script_path / config.kgbased_functional_network_rscript_path),
            "--nodes", to_r_path(config.outdir.parent.parent / "Node_info.tsv"),
            "--edges", to_r_path(config.outdir.parent.parent / "Edge_info.tsv"),
            "--module_info", to_r_path(config.outdir.parent / "DEP_community_detection_modules.tsv"),
            "--module_id", str(module_id),
            "--outdir", to_r_path(config.outdir.parent),
            "--enrich_fromType", getattr(config, 'identified_type', 'UNIPROT'),
            "--max_phosphoSite_displayed", str(max_phos_display),
            "--network_mode", network_mode,
            "--SEEDNODEID", seed_node,
            "--node_filtering", centrality_metric,
            "--network_layout", layout_algo,
            "--top_nodes_visualization_num", str(top_nodes),
            "--omics1_name", config.omics1_name,
            "--omics2_name", config.omics2_name,
            "--node_up", config.node_up,
            "--node_down", config.node_down,
            "--node_nonsig", config.node_nonsig,
            "--node_notdet", config.node_notdet,
            "--node_disease_border_color", config.node_disease_border_color,
            "--node_notdisease_border_color", config.node_notdisease_border_color,
            "--function_enrich_color_gradient_low", config.function_enrich_color_gradient_low,
            "--function_enrich_color_gradient_high", config.function_enrich_color_gradient_high,
            "--edge_type_colors", edge_colors_to_param(config.edge_type_colors),
            "--function_enrich_rscript_path", to_r_path(config.script_path / config.function_enrich_rscript_path)
        ]
        
        # Execute R visualization
        print(f"📋 R command for module {module_id}:\n{' '.join(r_cmd)}", file=sys.stderr)
        result = subprocess.run(r_cmd, capture_output=True, text=True)
        
        # Report execution status
        print(f"Module {module_id} - STDOUT:\n{result.stdout[:500]}...", file=sys.stderr)
        if result.stderr:
            print(f"Module {module_id} - STDERR:\n{result.stderr}", file=sys.stderr)
        
        if result.returncode == 0:
            print(f"✅ Module {module_id} visualization complete", file=sys.stderr)
        else:
            print(f"⚠️  Module {module_id} R script warnings (code {result.returncode})", file=sys.stderr)
    
    print(f"🎉 All visualizations saved to: {config.outdir}", file=sys.stderr)
    return config
