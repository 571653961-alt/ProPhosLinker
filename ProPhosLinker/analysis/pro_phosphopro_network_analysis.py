#!/usr/bin/env python3
"""
ProPhosLinker: Knowledge Graph-based Network Analysis Module (3.2)

This module performs comprehensive network analysis on protein-phosphoprotein 
functional interaction networks constructed from knowledge graphs. It provides:

Key Features:
- Network construction using NetworkX from node/edge TSV files
- Multiple network filtering modes (ALL, KS, LR, TF, DEP, DEPUP, DEPDOWN, SEEDNODE)
- Community detection using 6 algorithms (Louvain, Leiden, Infomap, FastGreedy, 
  Walktrap, Label Propagation)
- Integration with R visualization scripts
- Automated output directory management and result saving
- Support for phosphosite aggregation and disease association analysis

The pipeline processes proteomics/phosphoproteomics networks through:
1. Network construction and validation
2. Mode-specific subnetwork extraction
3. Community detection and modularity analysis
4. Integration with downstream R-based visualization

"""


import os
import sys
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
import igraph as ig
import leidenalg
from typing import Union, Dict, Literal
from pathlib import Path
import subprocess


class NetworkConfig:
    """
    Configuration class for protein-phosphoprotein network analysis.
    
    This class manages all parameters and data for network construction,
    analysis, and visualization. It handles input validation, data loading,
    and provides a centralized configuration for the entire pipeline.
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
        # Network data files (optional)
        nodefile: Union[str, Path, pd.DataFrame] = None,
        edgefile: Union[str, Path, pd.DataFrame] = None,
        # Module information (optional)
        module_info: Union[str, Path, pd.DataFrame] = None,
        # Omics names
        omics1_name: str = 'Pro',
        omics2_name: str = 'Phos',
        # Network construction parameters
        disease: str = '',
        kgbased_network_top_n: int = 6,
        # Network analysis parameters
        analysis_network_mode: Literal["ALL", "KS", "LR", "TF", "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'DEP',
        analysis_SEEDNODEID: str = None,
        comm_detection: Literal["louvain", "leiden", "infomap", "fastgreedy", "walktrap", "LPA"] = 'infomap',
        # Network visualization parameters
        top_nodes_visualization_num: int = 20,
        module_id: int = 0,
        max_phosphoSite_displayed: int = 5,
        visualization_network_mode: Literal["ALL", "KS", "LR", "TF", "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'ALL',
        visualization_SEEDNODEID: str = None,
        node_filtering: Literal["betweenness", "degree", "closeness", "harmonic", "eigenvector", "pagerank", "alpha", "hub", "authority"] = 'betweenness',
        network_layout: Literal["fr", "kk", "dh", "stress", "tree", "gem", "graphopt", "lgl", "circle", "grid"] = 'kk',
        oudir_default_name: str = "3.2KGbased_Functional_Network",
        # Color parameters
        node_up: str = "#B83A2D",
        node_down: str = "#657f68",
        node_nonsig: str = "#E3C79F",
        node_notdet: str = "grey90",
        node_disease_border_color: str = "#535c54",
        node_notdisease_border_color: str = "#acabab",
        function_enrich_color_gradient_low: str = "#175663",
        function_enrich_color_gradient_high: str = "#90362d",
        # Edge type color mapping
        edge_type_colors: Dict[str, str] = {
            "association": "#d7d6d6",           # Light gray
            "physical association": "#838181",  # Medium gray (same as binding)
            "binding": "#838181",               # Medium gray
            "direct interaction": "#d7d6d6",    # Light gray
            
            "activation": "#d62c0c",            # Orange-red
            "catalysis": "#7f00ff",             # Light purple
            "proximity": "#FF9301",             # Orange
            "reaction": "#114335",              # Dark green
            "phosphorylation": "#2FBE95",       # Blue-green
            "dephosphorylation": "#6B8E23",     # Olive green
            
            "ptmod": "#8C97D6",                 # Light purple
            "inhibition": "#0cb6d6",            # Dark blue
            "expression": "#FCF402",            # Yellow
            "regulation": "#e4dadd",            # Purple-gray
            
            "colocalization": "#4c95cd",        # Sky blue
            "covalent binding": "#716F74",      # Violet
            "ubiquitination": "#FF4500",        # Orange-red
            "multiRel": "#9a6728"               # Brown
        },
        # Script paths
        script_path: str = './scripts',
        kgbased_functional_network_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_visualization.R',
        kgbased_functional_network_community_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_community_visualisation.R',
        function_enrich_rscript_path: Union[str, Path] = '3.2functional_enrichment_function.R',
        
        # Neo4j database connection parameters
        uri: str = "bolt://localhost:7687",
        username: str = "neo4j",
        password: str = "neo4j"
    ):
        """Initialize NetworkConfig with provided parameters."""
        
        # Store all parameters as instance attributes
        self.profile = profile
        self.phosphoprofile = phosphoprofile
        self.prodiff = prodiff
        self.phosphoprodiff = phosphoprodiff
        self.phosphopro_pro_cor = phosphopro_pro_cor
        self.outdir = outdir
        self.nodefile = nodefile
        self.edgefile = edgefile
        self.module_info = module_info
        self.omics1_name = omics1_name
        self.omics2_name = omics2_name
        self.disease = disease
        self.kgbased_network_top_n = kgbased_network_top_n
        self.analysis_network_mode = analysis_network_mode
        self.analysis_SEEDNODEID = analysis_SEEDNODEID
        self.comm_detection = comm_detection
        self.top_nodes_visualization_num = top_nodes_visualization_num
        self.module_id = module_id
        self.max_phosphoSite_displayed = max_phosphoSite_displayed
        self.visualization_network_mode = visualization_network_mode
        self.visualization_SEEDNODEID = visualization_SEEDNODEID
        self.node_filtering = node_filtering
        self.network_layout = network_layout
        
        # Color parameters
        self.node_up = node_up
        self.node_down = node_down
        self.node_nonsig = node_nonsig
        self.node_notdet = node_notdet
        self.node_disease_border_color = node_disease_border_color
        self.node_notdisease_border_color = node_notdisease_border_color
        self.function_enrich_color_gradient_low = function_enrich_color_gradient_low
        self.function_enrich_color_gradient_high = function_enrich_color_gradient_high
        self.edge_type_colors = edge_type_colors
        
        # Script paths
        self.oudir_default_name = oudir_default_name
        self.script_path = script_path
        self.kgbased_functional_network_rscript_path = kgbased_functional_network_rscript_path
        self.kgbased_functional_network_community_rscript_path = kgbased_functional_network_community_rscript_path
        self.function_enrich_rscript_path = function_enrich_rscript_path
        
        # Neo4j connection parameters
        self.uri = uri
        self.username = username
        self.password = password

        # Validate parameters and load data
        self._validate_params()
        self._data_load()

    def _validate_params(self):
        """
        Validate input parameters and set up output directory.
        
        Raises:
            TypeError: If parameters have incorrect types
            ValueError: If output path is invalid
            FileNotFoundError: If input files don't exist
            RuntimeError: If directory creation fails
        """
        # Ensure outdir is a Path object
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        elif self.outdir is None:
            self.outdir = os.getcwd()
        elif not isinstance(self.outdir, Path):
            raise TypeError("outdir must be a string or Path object")
            
        # Check if output path exists and is writable
        if self.outdir.exists():
            if self.outdir.is_file():
                raise ValueError(f"Output path is a file, not a directory: {self.outdir}")
            if not os.access(self.outdir, os.W_OK):
                raise ValueError(f"Directory exists but is not writable: {self.outdir}")
                
        # Create output subdirectory with default name
        self.outdir = self.outdir / Path(self.oudir_default_name)
        output_dir = self.outdir
        
        # Create the output directory if it does not exist
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
            
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f"✅ Output directory ready: {output_dir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {output_dir}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")
        
        # Validate input files
        for attr in ['profile', 'phosphoprofile', 'prodiff', 'phosphoprodiff', 'phosphopro_pro_cor']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Path(value))
            elif not isinstance(value, (Path, pd.DataFrame)):
                raise TypeError(f"{attr} must be a string, Path object, or DataFrame")
                
            # Check file existence for Path objects
            if isinstance(value, Path) and not value.exists():
                raise FileNotFoundError(f"{value} does not exist")
            if isinstance(value, Path) and not value.is_file():
                raise ValueError(f"{value} is not a valid file")
                
            # Validate DataFrame content
            if isinstance(value, pd.DataFrame):
                if value.empty:
                    raise ValueError(f"{attr} DataFrame is empty")
                if not all(col in value.columns for col in ['ID', 'Name']):
                    raise ValueError(f"{attr} DataFrame must contain 'ID' and 'Name' columns")
        
        # Validate fold change threshold
        if self.FC <= 0:
            raise ValueError("FC (log2FC) value must be > 0")

    def _data_load(self):
        """
        Load and validate all input data files.
        
        Converts file paths to DataFrames and validates their content.
        Handles expression data, differential expression results,
        and phosphosite-protein mappings.
        """
        
        def datafile_load_validation(file):
            """
            Load and validate expression data files.
            
            Args:
                file: Path to file or DataFrame
                
            Returns:
                Validated DataFrame
                
            Raises:
                ValueError: If data is empty or has duplicate indices
            """
            if not isinstance(file, pd.DataFrame):
                # Read data from file
                df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
                df = df.apply(pd.to_numeric, errors='coerce')
            else:
                df = file
                
            if df.empty:
                raise ValueError(file.as_posix() + " Expression data is empty")
            if df.index.duplicated().any():
                raise ValueError(file.as_posix() + " Contains duplicate row indices")
            return df
        
        def difffile_load_validation(file):
            """
            Load and validate differential expression result files.
            
            Args:
                file: Path to file or DataFrame
                
            Returns:
                Validated DataFrame
            """
            if not isinstance(file, pd.DataFrame):
                df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
            else:
                df = file
            return df
        
        # Load all data files
        self.profile = datafile_load_validation(self.profile)
        self.phosphoprofile = datafile_load_validation(self.phosphoprofile)
        self.prodiff = difffile_load_validation(self.prodiff)
        self.phosphoprodiff = difffile_load_validation(self.phosphoprodiff)
        
        # Load phosphosite-protein mapping (as dictionary)
        if not isinstance(self.phosphopro_pro_cor, dict):
            self.phosphopro_pro_cor = pd.read_csv(
                self.phosphopro_pro_cor, sep='\t', header=None, 
                names=['value', 'key']
            ).set_index('key')['value'].to_dict()


def parameter_validation(nodefile, edgefile, output, analysis_network_mode):
    """
    Validate input parameters and prepare output directory.
    
    Args:
        nodefile: Path to node file
        edgefile: Path to edge file
        output: Output directory path
        analysis_network_mode: Network analysis mode
        
    Returns:
        Path: Validated output directory path
        
    Raises:
        FileExistsError: If input file doesn't exist
        ValueError: If input file is empty or output path is invalid
        RuntimeError: If directory creation fails
    """
    
    # Check if input file exists and is not empty
    def input_file_check(file_path):
        """Check if file exists and is not empty."""
        if not os.path.exists(file_path):
            raise FileExistsError(f"{file_path} Path does not exist")
        if os.path.getsize(file_path) == 0:
            raise ValueError(f"{file_path} File is empty!")
            
    # Prepare output directory
    def output_directory_check(output):
        """Create and validate output directory."""
        if output is None:
            base_path = os.getcwd()
        else:
            base_path = Path(output).expanduser().absolute()
            # Validate base path
            if base_path.exists():
                if base_path.is_file():
                    raise ValueError(f"Output path is a file, not a directory: {base_path}")
                if not os.access(base_path, os.W_OK):
                    raise ValueError(f"Directory exists but is not writable: {base_path}")
        
        # Create subdirectory for community detection results
        default_name = "3.3KGbased_Functional_Network/" + analysis_network_mode + "_Community_detection"
        output_dir = base_path / Path(default_name)
        
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f"✅ Output directory ready: {output_dir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {output_dir}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")
            
        return output_dir

    # Validate input files
    input_file_check(nodefile)
    input_file_check(edgefile)
    
    # Create and return output directory
    return output_directory_check(output)


def input_load_validation(nodefile, edgefile, network_mode, SEEDNODEID):
    """
    Load and validate input node and edge files.
    
    Args:
        nodefile: Path to node file
        edgefile: Path to edge file
        network_mode: Network filtering mode
        SEEDNODEID: Seed node ID (for SEEDNODE mode)
        
    Returns:
        tuple: (nodes_df, edges_df) DataFrames containing node and edge data
        
    Raises:
        ValueError: If data is empty, has duplicate indices, or seed node not found
    """
    
    def file_load_validation(file_path):
        """Load and validate a TSV file."""
        df = pd.read_csv(file_path, sep='\t')
        if df.empty:
            raise ValueError("Expression data is empty")
        if df.index.duplicated().any():
            raise ValueError("Contains duplicate row indices")
        return df
        
    # Load node and edge data
    nodes_df = file_load_validation(nodefile)
    edges_df = file_load_validation(edgefile)
    
    # Validate seed node if in SEEDNODE mode
    if network_mode == 'SEEDNODE' and SEEDNODEID:
        if SEEDNODEID not in nodes_df['NodeName'].tolist():
            raise ValueError(f"❌ {SEEDNODEID} Node does not exist in the network.")
            
    return nodes_df, edges_df


def network_construction(nodes_df, edges_df):
    """
    Construct a NetworkX graph from node and edge data.
    
    This function:
        1. Aggregates multiple phosphosite records per protein
        2. Creates nodes with all attributes
        3. Adds edges with interaction types
        
    Args:
        nodes_df: DataFrame with node information
        edges_df: DataFrame with edge information
        
    Returns:
        nx.Graph: Constructed network graph
    """
    
    # Preprocess nodes by aggregating multiple phosphosite records
    def preprocess_nodes(df):
        """
        Group multiple rows for the same NodeName and aggregate phosphosite information.
        
        Args:
            df: Raw node DataFrame
            
        Returns:
            DataFrame: Processed nodes with aggregated phosphosite data
        """
        
        # Use defaultdict to collect data for each node
        grouped = defaultdict(lambda: {
            # Protein information
            'Pro_type': None,
            'Pro_FC': None,
            'Pro_FDR': None,
            'Pro_class': None,
            'DiseaseRelated': None,
            'phosphoSite_detected_num': 0,
            # Phosphorylation information
            'phosphoSites': [],  # Store all phosphosite information
        })
        
        # Define column names
        protein_cols = ['Pro_type', 'Pro_FC', 'Pro_FDR', 'Pro_class', 'DiseaseRelated', 'phosphoSite_detected_num']
        phospho_cols = ['phosphoSite_index', 'phosphoprotein_ID', 'Phosphopro_type', 'Phosphopro_FC', 'Phosphopro_FDR', 'Phosphopro_class']

        # Collect all attributes
        for _, row in df.iterrows():
            node_name = row['NodeName']
            entry = grouped[node_name]
            
            # Update protein information (take first non-null value)
            for col in protein_cols:
                if col in row and (entry[col] is None or pd.isna(entry[col])) and pd.notna(row[col]):
                    entry[col] = row[col]
        
            # Collect phosphosite information (if available)
            phospho_data = {}
            for col in phospho_cols:
                if col in row and pd.notna(row[col]):
                    phospho_data[col] = row[col]
            if phospho_data:
                entry['phosphoSites'].append(phospho_data)

        # Convert to DataFrame
        processed_data = []
        for node_name, data in grouped.items():
            # Create node data
            node_data = {
                'NodeName': node_name,
                'num_phosphoSites': len(data['phosphoSites'])
            }
            # Add protein attributes
            for col in protein_cols:
                node_data[col] = data[col]
            # Add phosphosite information
            node_data['phosphoSites'] = data['phosphoSites']
            processed_data.append(node_data)

        return pd.DataFrame(processed_data)
    
    # Process node data
    processed_nodes = preprocess_nodes(nodes_df)
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes with attributes
    for _, row in processed_nodes.iterrows():
        node_id = row['NodeName']
        attrs = row.to_dict()
        # Convert numpy types to Python native types
        for k, v in attrs.items():
            if isinstance(v, np.generic):
                attrs[k] = v.item() if np.isscalar(v) else v.tolist()
        G.add_node(node_id, **attrs)
    
    # Add edges with attributes
    if not edges_df.empty:
        for _, row in edges_df.iterrows():
            source = row['Source']
            target = row['Target']
            # Ensure both nodes exist in the graph
            if source in G and target in G:
                edges_attrs = {'Interaction_type': row['Interaction_Type']}
                G.add_edge(source, target, **edges_attrs)

    return G


def network_mode_extraction(G_raw, nodes_df, network_mode, SEEDNODEID):
    """
    Extract subnetwork based on specified filtering mode.
    
    Supported modes:
        - ALL: Keep entire network
        - KS: Keep only kinase-substrate interactions (edges containing 'phospho')
        - LR: Keep only ligand-receptor interactions (edges containing 'binding')
        - TF: Keep only transcription factor interactions (activation, regulation, etc.)
        - DEP: Keep nodes with differential expression (UP or DOWN)
        - DEPUP: Keep only upregulated nodes
        - DEPDOWN: Keep only downregulated nodes
        - SEEDNODE: Extract subnetwork around seed node(s)
    
    Args:
        G_raw: Original network graph
        nodes_df: DataFrame with node information
        network_mode: Filtering mode
        SEEDNODEID: Seed node ID (for SEEDNODE mode)
        
    Returns:
        nx.Graph: Filtered subnetwork
    """
    
    if network_mode == 'ALL':
        # Keep entire network
        G = G_raw.copy()
        
    elif network_mode == 'KS':
        # Filter to keep only kinase-substrate interactions
        G = G_raw.edge_subgraph([(u, v) for u, v, d in G_raw.edges(data=True)
                                 if isinstance(d.get('Interaction_type'), str) and 'phospho' in d.get('Interaction_type')])
        
    elif network_mode == 'LR':
        # Filter to keep only ligand-receptor interactions
        G = G_raw.edge_subgraph([(u, v) for u, v, d in G_raw.edges(data=True)
                                 if isinstance(d.get('Interaction_type'), str) and 'binding' in d.get('Interaction_type')])
        
    elif network_mode == 'TF':
        # Filter to keep only transcription factor related interactions
        keywords = ['activation', 'regulation', 'repression', 'inhibition']
        G = G_raw.edge_subgraph([(u, v) for u, v, d in G_raw.edges(data=True)
                                 if isinstance(d.get('Interaction_type'), str) and any(kw in d.get('Interaction_type') for kw in keywords)])
        
    elif network_mode == 'DEP':
        # Keep differentially expressed nodes (UP or DOWN in either protein or phosphosite)
        deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('UP|DOWN', na=False)) |
                             (nodes_df['Phosphopro_class'].str.contains('UP|DOWN', na=False))]['NodeName'].unique().tolist()
        G = G_raw.subgraph(deg_nodes).copy()
        
    elif network_mode == 'DEPUP':
        # Keep only upregulated nodes
        deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('UP', na=False)) |
                             (nodes_df['Phosphopro_class'].str.contains('UP', na=False))]['NodeName'].unique().tolist()
        G = G_raw.subgraph(deg_nodes).copy()
        
    elif network_mode == 'DEPDOWN':
        # Keep only downregulated nodes
        deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('DOWN', na=False)) |
                             (nodes_df['Phosphopro_class'].str.contains('DOWN', na=False))]['NodeName'].unique().tolist()
        G = G_raw.subgraph(deg_nodes).copy()
        
    elif network_mode == 'SEEDNODE' and SEEDNODEID:
        # Extract subnetwork around seed node (include first and second-degree neighbors)
        SEEDNODE_neighbors = list(G_raw.neighbors(SEEDNODEID))
        
        if len(SEEDNODE_neighbors) < 10:
            # Include second-degree neighbors if first-degree are insufficient
            second_neighbors = set()
            for node in SEEDNODE_neighbors:
                neighbors = set(G_raw.neighbors(node))
                second_neighbors.update(neighbors)
            seed_neighbors_subnodes = list(set(SEEDNODE_neighbors + list(second_neighbors) + [SEEDNODEID]))
        else:
            seed_neighbors_subnodes = list(set(SEEDNODE_neighbors + [SEEDNODEID]))
            
        G = G_raw.subgraph(seed_neighbors_subnodes).copy()
        
    else:
        # Default: keep entire network
        G = G_raw.copy()
        
    return G


def network_community_detection(G, network_mode, comm_detection_type, out_dir):
    """
    Perform community detection on the network using specified algorithm.
    
    Supported algorithms:
        - louvain: Louvain community detection
        - leiden: Leiden algorithm (higher quality communities)
        - infomap: Infomap algorithm (flow-based)
        - fastgreedy: Fast greedy modularity optimization
        - walktrap: Walktrap algorithm (random walks)
        - LPA: Label Propagation Algorithm
    
    Args:
        G: NetworkX graph
        network_mode: Network filtering mode
        comm_detection_type: Community detection algorithm to use
        out_dir: Output directory for saving results
        
    Returns:
        DataFrame: Module information with ModuleID, ModuleSize, and ModuleItems
    """
    
    output_path = Path(out_dir) / f"{network_mode}_community_detection_modules.tsv"
    
    # Prepare list to store module data
    module_data = []
    
    # Identify and remove isolated nodes
    isolated_nodes = [node for node, degree in dict(G.degree()).items() if degree == 0]
    
    # Create subgraph without isolated nodes
    G_main = G.copy()
    G_main.remove_nodes_from(isolated_nodes)
    
    # Perform community detection on the main graph (not for SEEDNODE mode)
    if network_mode != 'SEEDNODE':
        if G_main.number_of_edges() > 0:
            # Convert to igraph for community detection
            G_ig = ig.Graph.from_networkx(G_main)
            
            # Apply selected algorithm
            if comm_detection_type == 'louvain':
                partition = G_ig.community_multilevel()
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            elif comm_detection_type == 'leiden':
                partition = leidenalg.find_partition(G_ig, leidenalg.RBConfigurationVertexPartition, resolution_parameter=1.0)
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            elif comm_detection_type == 'infomap':
                partition = G_ig.community_infomap()
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            elif comm_detection_type == 'fastgreedy':
                dendrogram = G_ig.community_fastgreedy()
                partition = dendrogram.as_clustering()
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            elif comm_detection_type == 'walktrap':
                dendrogram = G_ig.community_walktrap(steps=4)  # steps parameter controls random walk length
                partition = dendrogram.as_clustering()  # Get optimal partition (max modularity)
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            elif comm_detection_type == 'LPA':
                partition = G_ig.community_label_propagation()
                nx.set_node_attributes(G_main, partition.membership, 'community')
                
            else:
                # No community detection
                nx.set_node_attributes(G_main, np.zeros(len(G_main.nodes)), 'community')
                partition = None
                
            # modularity = G_ig.modularity(partition.membership)  # Optionally calculate modularity
        else:
            partition = None
            # modularity = 0.0

        # Analyze community structure (main graph)
        community_dict = {}
        
        # Ensure partition exists
        if partition:
            # partition.membership is a list ordered by G_ig.vs indices
            for node, comm_id in zip(G_main.nodes, partition.membership):
                if comm_id not in community_dict:
                    community_dict[comm_id] = []
                community_dict[comm_id].append(node)
        else:
            # If no edges, all nodes are isolated
            for node in G_main.nodes:
                if 0 not in community_dict:
                    community_dict[0] = []
                community_dict[0].append(node)
                
        # Sort communities by size (descending)
        sorted_communities = sorted(community_dict.items(), key=lambda x: len(x[1]), reverse=True)

        # Collect community information
        for comm_id, members in sorted_communities:
            module_data.append({
                'ModuleID': comm_id,
                'ModuleSize': len(members),
                'ModuleItems': ';'.join(members)
            })
            
        # Handle isolated nodes
        if not partition:  # No partition information
            module_data.append({
                'ModuleID': -1,
                'ModuleSize': len(isolated_nodes),
                'ModuleItems': ';'.join(isolated_nodes)
            })
        elif -1 not in community_dict:  # Partition exists but doesn't include isolated nodes
            module_data.append({
                'ModuleID': -1,
                'ModuleSize': len(isolated_nodes),
                'ModuleItems': ';'.join(isolated_nodes)
            })

        # Create DataFrame and save
        f_Module_info = pd.DataFrame(module_data, columns=['ModuleID', 'ModuleSize', 'ModuleItems'])
        f_Module_info.to_csv(output_path, sep='\t', index=False)

    else:
        # For SEEDNODE mode, skip community detection and just output all nodes
        output_path = Path(out_dir) / f"{network_mode}_community_detection_modules.tsv"
        module_data.append({
            'ModuleID': 0,
            'ModuleSize': len(G_main.nodes),
            'ModuleItems': ';'.join(G_main.nodes)
        })
        
        # Create DataFrame and save
        f_Module_info = pd.DataFrame(module_data, columns=['ModuleID', 'ModuleSize', 'ModuleItems'])
        f_Module_info.to_csv(output_path, sep='\t', index=False)
            
    return f_Module_info


def NetworkAnalysis(config: NetworkConfig):
    """
    Main network analysis pipeline.
    
    This function orchestrates the entire network analysis workflow:
        1. Sets up output directory
        2. Constructs the network
        3. Applies mode-based filtering
        4. Performs community detection
        5. Calls R scripts for visualization
    
    Args:
        config: NetworkConfig object with all parameters
    """
    
    # Set up output directory for this analysis mode
    config.outdir = config.outdir / Path(config.analysis_network_mode + '_community_detection')
    if not config.outdir.exists():
        config.outdir.mkdir(parents=True, exist_ok=True)
        
    # Step 1: Network construction
    G_raw = network_construction(config.nodefile, config.edgefile)
    
    # Step 2: Mode-based extraction
    G = network_mode_extraction(G_raw, config.nodefile, config.analysis_network_mode, config.analysis_SEEDNODEID)
    
    # Step 3: Community detection
    config.module_info = network_community_detection(G, config.analysis_network_mode, config.comm_detection, config.outdir)
    
    # Step 4: Visualization using R script
    network_mode = config.visualization_network_mode

    def to_r_path(path):
        """Convert Windows path to R-compatible path with forward slashes."""
        return str(path).replace("\\", "/")
        
    curr_outdir = config.outdir.parent
    
    # Prepare R script command
    cmd = [
        "Rscript",
        to_r_path(str(config.script_path / Path(config.kgbased_functional_network_community_rscript_path))),
        "--edge_file", to_r_path(curr_outdir / Path('Edge_diff_info.tsv')),
        "--community_info", to_r_path(curr_outdir / Path('DEP_community_detection') / Path('DEP_community_detection_modules.tsv')),
        "--outdir", to_r_path(curr_outdir / Path('DEP_community_detection')),
        "--top_module_num", str(config.kgbased_network_top_n)
    ]

    full_cmd = cmd
    # Execute R script
    print(f"📋 Executing command:\n{' '.join(full_cmd)}")
    result = subprocess.run(full_cmd, capture_output=True, text=True)

    # Print output for debugging
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)



# import argparse
# import os
# import sys
# from pathlib import Path
# import networkx as nx
# import numpy as np
# import pandas as pd
# from collections import defaultdict
# import igraph as ig
# import leidenalg
# from typing import Union, Dict,Literal
# from pathlib import Path
# import subprocess


# class NetworkConfig:
#     def __init__(
#         self,
#         *,
#         #input file 
#         profile : Union[str, Path, pd.DataFrame],
#         phosphoprofile : Union[str, Path, pd.DataFrame],
#         prodiff : Union[str, Path, pd.DataFrame],
#         phosphoprodiff : Union[str, Path, pd.DataFrame],
#         phosphopro_pro_cor : Union[str, Path, Dict],
#         # output directory
#         outdir : Union[str, Path],
#         # the parameter of the network data
#         nodefile : Union[str, Path, pd.DataFrame] = None,
#         edgefile : Union[str, Path, pd.DataFrame] = None,
#         # 
#         module_info : Union[str, Path, pd.DataFrame] = None,
#         #omics name
#         omics1_name : str = 'Pro',
#         omics2_name : str = 'Phos',
#         # the parameter of the network constructor 
#         disease : str = '',
#         kgbased_network_top_n : int = 6,
#         # the parameter of the network analysis
#         analysis_network_mode : Literal["ALL",    # 3 catalogs
#                         "KS","LR","TF",
#                         "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'DEP',
#         analysis_SEEDNODEID : str = None,
#         comm_detection : Literal["louvain","leiden","infomap","fastgreedy","walktrap", "LPA"] ='infomap',
#         #  the parameter of the network visualization
#         top_nodes_visualization_num: int = 20,
#         module_id: int = 0,
#         max_phosphoSite_displayed: int = 5,
#         visualization_network_mode : Literal["ALL",    # 3 catalogs
#                         "KS","LR","TF",
#                         "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'ALL',
#         visualization_SEEDNODEID : str = None,
#         node_filtering : Literal["betweenness", "degree",  "closeness", "harmonic","eigenvector", "pagerank", "alpha", "hub",  "authority"] = 'betweenness',
#         network_layout : Literal["fr","kk","dh","stress","tree","gem","graphopt", "lgl","circle","grid"] = 'kk',
#         oudir_default_name : str = "3.2KGbased_Functional_Network",
#         node_up :str = "#B83A2D",
#         node_down :str = "#657f68",
#         node_nonsig :str = "#E3C79F",
#         node_notdet :str = "grey90",
#         node_disease_border_color :str = "#535c54",
#         node_notdisease_border_color :str = "#acabab",
#         function_enrich_color_gradient_low :str = "#175663",
#         function_enrich_color_gradient_high :str = "#90362d",
#         # 边类型颜色映射
#         edge_type_colors: Dict[str, str] = {
#             "association": "#d7d6d6",           # 浅灰
#             "physical association": "#838181",  # 浅灰（与 binding 相同）
#             "binding": "#838181",               # 浅灰
#             "direct interaction": "#d7d6d6",    # 浅灰
            
#             "activation": "#d62c0c",            # 橙红色
#             "catalysis": "#7f00ff",             # 淡紫色
#             "proximity": "#FF9301",             # 橙色
#             "reaction": "#114335",              # 深绿色
#             "phosphorylation": "#2FBE95",       # 蓝绿色
#             "dephosphorylation": "#6B8E23",     # 橄榄绿
            
#             "ptmod": "#8C97D6",                 # 淡紫色
#             "inhibition": "#0cb6d6",            # 深蓝色
#             "expression": "#FCF402",            # 黄色
#             "regulation": "#e4dadd",            # 紫灰色
            
#             "colocalization": "#4c95cd",        # 天蓝色
#             "covalent binding": "#716F74",      # 紫罗兰色
#             "ubiquitination": "#FF4500",        # 橙红色
#             "multiRel": "#9a6728"               # 棕色
#         },
#         #
#         script_path : str = './scripts',
#         kgbased_functional_network_rscript_path : Union[str, Path] = '3.2functionnal_interaction_network_visualization.R',
#         kgbased_functional_network_community_rscript_path : Union[str, Path] = '3.2functionnal_interaction_network_community_visualisation.R',
#         function_enrich_rscript_path : Union[str, Path] = '3.2functional_enrichment_function.R',
        
#         ####
#         uri : str =  "bolt://localhost:7687",
#         username : str =  "neo4j",
#         password : str =  "neo4j"
#         ):
#         self.profile = profile
#         self.phosphoprofile = phosphoprofile
#         self.prodiff = prodiff
#         self.phosphoprodiff = phosphoprodiff
#         self.phosphopro_pro_cor = phosphopro_pro_cor
#         self.outdir = outdir
#         self.nodefile = nodefile
#         self.edgefile = edgefile
#         self.module_info = module_info
#         self.omics1_name = omics1_name
#         self.omics2_name = omics2_name
#         self.FC = FC
#         self.disease = disease
#         self.kgbased_network_top_n = kgbased_network_top_n
#         self.analysis_network_mode = analysis_network_mode
#         self.analysis_SEEDNODEID = analysis_SEEDNODEID
#         self.comm_detection = comm_detection
#         self.top_nodes_visualization_num = top_nodes_visualization_num
#         self.module_id = module_id
#         self.max_phosphoSite_displayed = max_phosphoSite_displayed
#         self.visualization_network_mode = visualization_network_mode
#         self.visualization_SEEDNODEID = visualization_SEEDNODEID
#         self.node_filtering = node_filtering
#         self.network_layout = network_layout
#         #
#         self.node_up = node_up
#         self.node_down = node_down
#         self.node_nonsig = node_nonsig
#         self.node_notdet = node_notdet
#         self.node_disease_border_color = node_disease_border_color
#         self.node_notdisease_border_color = node_notdisease_border_color
#         self.function_enrich_color_gradient_low = function_enrich_color_gradient_low
#         self.function_enrich_color_gradient_high = function_enrich_color_gradient_high
#         # 边类型颜色映射
#         self.edge_type_colors = edge_type_colors
#         #
#         self.oudir_default_name = oudir_default_name
#         self.script_path = script_path
#         self.kgbased_functional_network_rscript_path  = kgbased_functional_network_rscript_path
#         self.kgbased_functional_network_community_rscript_path  = kgbased_functional_network_community_rscript_path
#         self.function_enrich_rscript_path = function_enrich_rscript_path
#         #
#         self.uri = uri
#         self.username = username
#         self.password = password

#         self._validate_params()
#         self._data_load()

#     def _validate_params(self):
#         # Ensure outdir is a Path object
#         if isinstance(self.outdir, str):
#             self.outdir = Path(self.outdir)
#         elif self.outdir is None:
#             self.outdir = os.getcwd()
#         elif not isinstance(self.outdir, Path):
#             raise TypeError("outdir must be a string or Path object")
#         if self.outdir.exists():
#             if self.outdir.is_file():
#                 raise ValueError(f"指定输出路径时文件而非目录：{self.outdir}")
#             if not os.access(self.outdir, os.W_OK):
#                 raise ValueError(f"目录已存在且不允许覆盖: {self.outdir}")
#         self.outdir = self.outdir/ Path(self.oudir_default_name)
#         output_dir = self.outdir
#         # Create the output directory if it does not exist
#         if not self.outdir.exists():
#             self.outdir.mkdir(parents=True, exist_ok=True)
#         try:
#             output_dir.mkdir(parents=True,exist_ok=True)
#             print(f"✅ 输出目录就绪: {output_dir}", file=sys.stderr)
#         except PermissionError:
#              RuntimeError(f"无权限创建目录: {output_dir}")
#         except Exception as e:
#             raise RuntimeError(f"目录创建失败: {str(e)}")
        
#         # Ensure all input files are Path objects or DataFrames
#         for attr in ['profile', 'phosphoprofile', 'prodiff', 'phosphoprodiff', 'phosphopro_pro_cor']:
#             value = getattr(self, attr)
#             if isinstance(value, str):
#                 setattr(self, attr, Path(value))
#             elif not isinstance(value, (Path, pd.DataFrame)):
#                 raise TypeError(f"{attr} must be a string, Path object, or DataFrame")
#             # Ensure the input files exist if they are Path objects
#             if isinstance(value, Path) and not value.exists():
#                 raise FileNotFoundError(f"{value} does not exist")
#             # Ensure the input files are readable if they are Path objects
#             if isinstance(value, Path) and not value.is_file():
#                 raise ValueError(f"{value} is not a valid file")
#             # Ensure the input files are DataFrames if they are not Path objects
#             if isinstance(value, pd.DataFrame):
#                 if value.empty:
#                     raise ValueError(f"{attr} DataFrame is empty")
#                 if not all(col in value.columns for col in ['ID', 'Name']):
#                     raise ValueError(f"{attr} DataFrame must contain 'ID' and 'Name' columns")
    
#         if self.FC<= 0 :
#             raise ValueError("FC(log2FC)值 必须在 >0")
        

#     def _data_load(self):
#         def datafile_load_validation(file):
#             if not isinstance(file, pd.DataFrame):
#                     # 读取数据
#                     df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
#                     df = df.apply(pd.to_numeric, errors='coerce')
#             else:
#                 df = file
#             if df.empty:
#                 raise ValueError(file.as_posix()+"表达谱数据为空")
#             if df.index.duplicated().any():
#                 raise ValueError(file.as_posix()+"存在重复的行索引")
#             return df
        
#         def difffile_load_validation(file):
#             if not isinstance(file, pd.DataFrame):
#                     df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
#             else:
#                 df = file
#             return df
        
#         self.profile = datafile_load_validation(self.profile)
#         self.phosphoprofile = datafile_load_validation(self.phosphoprofile)
#         self.prodiff = difffile_load_validation(self.prodiff)
#         self.phosphoprodiff = difffile_load_validation(self.phosphoprodiff)
#         if not isinstance(self.phosphopro_pro_cor, dict):
#             self.phosphopro_pro_cor = pd.read_csv(self.phosphopro_pro_cor,sep='\t',header=None,names=['value','key']).set_index('key')['value'].to_dict()




# def parameter_validation(nodefile,edgefile,output,analysis_network_mode):
#     # 输入文件是否存在、输入文件大小是否为零、格式是否为.tsv
#     def input_file_check(file_path):
#         if not os.path.exists(file_path):
#             raise FileExistsError(f"{file_path} 该路径不存在")
#         if os.path.getsize(file_path) == 0:
#             raise ValueError(f"{file_path} 该文件内容为空！")
#     def output_directory_check(output):
#         if output is None:
#             base_path = os.getcwd()
#         else:
#             base_path = Path(output).expanduser().absolute()
#             # 验证基础路径类型
#             if base_path.exists():
#                 if base_path.is_file():
#                     raise ValueError(f"指定输出路径时文件而非目录：{base_path}")
#                 if not os.access(base_path, os.W_OK):
#                     raise ValueError(f"目录已存在且不允许覆盖: {base_path}")
#         # 验证输出文件没问题，创建results/Network_analysis子目录

#         default_name = "3.3KGbased_Functional_Network/"+analysis_network_mode+"_Community_detection"
#         output_dir = base_path / Path(default_name)
#         try:
#             output_dir.mkdir(parents=True,exist_ok=True)
#             print(f"✅ 输出目录就绪: {output_dir}", file=sys.stderr)
#         except PermissionError:
#             raise RuntimeError(f"无权限创建目录: {output_dir}")
#         except Exception as e:
#             raise RuntimeError(f"目录创建失败: {str(e)}")
#         return output_dir

#     input_file_check(nodefile)
#     input_file_check(edgefile)
#     #输入文件夹路径是否存在，没有就创建
#     return output_directory_check(output)

# def input_load_validation(nodefile,edgefile,network_mode,SEEDNODEID):
#     def file_load_validation(file_path):
#         df = pd.read_csv(file_path, sep='\t')
#         if df.empty:
#             raise ValueError("表达谱数据为空")
#         if df.index.duplicated().any():
#             raise ValueError("存在重复的行索引")
#         return df
#     nodes_df = file_load_validation(nodefile)
#     edges_df = file_load_validation(edgefile)
#     if network_mode == 'SEEDNODE' and SEEDNODEID:
#         if SEEDNODEID not in nodes_df['NodeName'].tolist():
#             raise ValueError(f"❌   {SEEDNODEID} 该节点并不在存在于网络中。")
#     return nodes_df,edges_df

# def network_construction(nodes_df,edges_df):
#     # 合并同一NodeName的多行记录，聚合磷酸化信息
#     def preprocess_nodes(df):
#         # 分组收集数据
#         grouped = defaultdict(lambda: {
#             # 蛋白质信息
#         'Pro_type': None,
#         'Pro_FC': None,
#         'Pro_FDR': None,
#         'Pro_class': None,
#         'DiseaseRelated': None,
#         'phosphoSite_detected_num': 0,
#         # 磷酸化信息
#         'phosphoSites': [],  # 存储所有磷酸化位点信息
#     })
#         # 列名定义
#         protein_cols = ['Pro_type', 'Pro_FC', 'Pro_FDR', 'Pro_class', 'DiseaseRelated', 'phosphoSite_detected_num']
#         phospho_cols = ['phosphoSite_index', 'phosphoprotein_ID', 'Phosphopro_type', 'Phosphopro_FC', 'Phosphopro_FDR', 'Phosphopro_class']

#         # 收集所有属性
#         for _,row in df.iterrows():
#             node_name = row['NodeName']
#             entry = grouped[node_name]
#             # 更新蛋白质信息（取第一个非空值）
#             for col in protein_cols:
#                 if col in row and (entry[col] is None or pd.isna(entry[col])) and pd.notna(row[col]):
#                     entry[col] = row[col]
        
#             # 收集磷酸化位点信息（如果有）
#             phospho_data = {}
#             for col in phospho_cols:
#                 if col in row and pd.notna(row[col]):
#                     phospho_data[col] = row[col]
#             if phospho_data:
#                 entry['phosphoSites'].append(phospho_data)

#         # 转换为 DataFrame
#         processed_data = []
#         for node_name,data in grouped.items():
#             #创建节点数据
#             node_data = {
#                 'NodeName' : node_name,
#                 'num_phosphoSites' :len(data['phosphoSites'])
#             }
#             #添加蛋白属性
#             for col in protein_cols:
#                 node_data[col] = data[col]
#             # 添加磷酸化信息
#             node_data['phosphoSites'] = data['phosphoSites']
#             processed_data.append(node_data)

#         return pd.DataFrame(processed_data)
#     # 处理数据节点
#     processed_nodes = preprocess_nodes(nodes_df)
#     # 创建图 
#     G = nx.Graph()
#     #添加节点
#     for _,row in processed_nodes.iterrows():
#         node_id = row['NodeName']
#         attrs = row.to_dict()
#         #转换为可能的numpy类型，为python 原生类型
#         for k,v in attrs.items():
#             if isinstance(v,np.generic):
#                 attrs[k] = v.item() if np.isscalar(v) else v.tolist()
#         G.add_node(node_id,**attrs)
#     #添加边
#     if not edges_df.empty:
#         for _,row in edges_df.iterrows():
#             source = row['Source']
#             target = row['Target']
#             # 确保边两端的节点存在
#             if source in G and target in G:
#                 edges_attrs = {'Interaction_type': row['Interaction_Type']}
#             G.add_edge(source,target,**edges_attrs)

#     return G

# def network_mode_extraction(G_raw,nodes_df,network_mode,SEEDNODEID):
#     if network_mode == 'ALL':
#         G = G_raw.copy()
#     elif network_mode == 'KS':
#         #筛选出G中 'kinase-substrate' 的子图，边包含'phospho'
#         G = G_raw.edge_subgraph([(u,v) for u,v,d in G_raw.edges(data=True)
#                                  if isinstance(d.get('Interaction_type'),str) and 'phospho' in d.get('Interaction_type')])
#     elif network_mode == 'LR':
#         #筛选出G中 'ligand-receptor' 的子图，边包含'binding'
#         G = G_raw.edge_subgraph([(u,v) for u,v,d in G_raw.edges(data =True)
#                                  if isinstance(d.get('Interaction_type'),str) and 'binding' in d.get('Interaction_type')])
#     elif network_mode == 'TF':
#         keywords = ['activation', 'regulation', 'repression', 'inhibition']
#         G = G_raw.edge_subgraph([(u,v) for u,v,d in G_raw.edges(data=True)
#                                  if isinstance(d.get('Interaction_type',str) and any(kw in d.get('Interaction_type') for kw in keywords))])
#     elif network_mode == 'DEP':
#         deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('UP|DOWN', na=False)) |(nodes_df['Phosphopro_class'].str.contains('UP|DOWN', na=False))]['NodeName'].unique().tolist()
#         G = G_raw.subgraph(deg_nodes).copy()
#     elif network_mode == 'DEPUP':
#         deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('UP', na=False)) |(nodes_df['Phosphopro_class'].str.contains('UP', na=False))]['NodeName'].unique().tolist()
#         G = G_raw.subgraph(deg_nodes).copy()
#     elif network_mode == 'DEPDOWN':
#         deg_nodes = nodes_df[(nodes_df['Pro_class'].str.contains('DOWN', na=False)) |(nodes_df['Phosphopro_class'].str.contains('DOWN', na=False))]['NodeName'].unique().tolist()
#         G = G_raw.subgraph(deg_nodes).copy()
#     elif network_mode == 'SEEDNODE' and SEEDNODEID:
#         SEEDNODE_neighbors = list(G_raw.neighbors(SEEDNODEID))
#         if(len(SEEDNODE_neighbors) < 10):
#             second_neighbors = set()
#             for node in SEEDNODE_neighbors:
#                 neighbors = set(G_raw.neighbors(node))
#                 second_neighbors.update(neighbors)
#             seed_neighbors_subnodes = list(set(SEEDNODE_neighbors+list(second_neighbors)+[SEEDNODEID]))
#         else:
#             seed_neighbors_subnodes = list(set(SEEDNODE_neighbors+[SEEDNODEID]))
#         G = G_raw.subgraph(seed_neighbors_subnodes).copy()
        
#     else:
#         G = G_raw.copy()
#     return G

# def network_community_detection(G,network_mode,comm_detection_type,out_dir):
#     output_path = Path(out_dir) / f"{network_mode}_community_detection_modules.tsv"
#     # 准备存储社区数据的列表
#     module_data = []
#     # 识别并移除孤立节点
#     isolated_nodes = [node for node,degree in dict(G.degree()).items() if degree == 0]
#     # 创建子图，移除孤立节点
#     G_main = G.copy()
#     G_main.remove_nodes_from(isolated_nodes)
#     # 只在右边的图上进行社群检测
#     if network_mode != 'SEEDNODE':
#         if G_main.number_of_edges() > 0:
#             G_ig = ig.Graph.from_networkx(G_main)
#             if comm_detection_type == 'louvain':
#                 partition = G_ig.community_multilevel()
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             elif comm_detection_type == 'leiden':
#                 partition = leidenalg.find_partition(G_ig, leidenalg.RBConfigurationVertexPartition, resolution_parameter=1.0)
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             elif comm_detection_type == 'infomap':
#                 partition = G_ig.community_infomap()
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             elif comm_detection_type == 'fastgreedy':
#                 dendrogram = G_ig.community_fastgreedy()
#                 partition = dendrogram.as_clustering()
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             elif comm_detection_type == 'walktrap':
#                 dendrogram = G_ig.community_walktrap(steps=4)  # steps参数控制随机游走步数
#                 # 获取最优划分（模块度最大的划分）
#                 partition = dendrogram.as_clustering()
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             elif comm_detection_type == 'LPA':
#                 partition = G_ig.community_label_propagation()
#                 nx.set_node_attributes(G_main,partition.membership,'community')
#             # elif comm_detection_type == 'spinglass':
#             #     partition = G_ig.community_spinglass()
#             else:
#                 nx.set_node_attributes(G_main,np.zeros(len(G_main.nodes)),'community')
#                 partition = None
#             # modularity = G_ig.modularity(partition.membership)
#         else:
#             partition = None
#             # modularity = 0.0

#         # 分析社群结构 （主图）
#         community_dict = {}
#         #确保划分结果存在
#         if partition:
#             # partition.membership 是一个列表，顺序与G_ig.vs 索引一致
#             for node,comm_id in zip(G_main.nodes,partition.membership):
#                 if comm_id not in community_dict:
#                     community_dict[comm_id] = []
#                 community_dict[comm_id].append(node)
#         else:
#             # 如果没有边，所有节点都是孤立的
#             for node in G_main.nodes:
#                 if 0 not in community_dict:
#                     community_dict[0] = []
#                 community_dict[0].append(node)
#         # 按社群成员数量排序
#         sorted_communities = sorted(community_dict.items(),key=lambda x:len(x[1]),reverse=True)

#         # 收集社区信息
#         for comm_id, members in sorted_communities:
#             module_data.append({
#                 'ModuleID': comm_id,
#                 'ModuleSize': len(members),
#                 'ModuleItems': ';'.join(members)
#             })
#         # 处理孤立节点
#         if not partition:  # 如果没有分区信息
#             module_data.append({
#                 'ModuleID': -1,
#                 'ModuleSize': len(isolated_nodes),
#                 'ModuleItems': ';'.join(isolated_nodes)
#             })
#         elif -1 not in community_dict:  # 如果分区存在但未包含孤立节点
#             module_data.append({
#                 'ModuleID': -1,
#                 'ModuleSize': len(isolated_nodes),
#                 'ModuleItems': ';'.join(isolated_nodes)
#             })

#         # 创建DataFrame
#         f_Module_info = pd.DataFrame(module_data, columns=['ModuleID', 'ModuleSize', 'ModuleItems'])
#         f_Module_info.to_csv(output_path, sep='\t', index=False)

#     else:
#         #不做社群检测，直接写入所有节点信息 
#         output_path = Path(out_dir) / f"{network_mode}_community_detection_modules.tsv"
#         module_data.append({
#                 'ModuleID': 0,
#                 'ModuleSize': len(G_main.nodes),
#                 'ModuleItems': ';'.join(G_main.nodes)
#             })
#         # 创建DataFrame
#         f_Module_info = pd.DataFrame(module_data, columns=['ModuleID', 'ModuleSize', 'ModuleItems'])
#         f_Module_info.to_csv(output_path, sep='\t', index=False)
            
#     return f_Module_info
            



# def NetworkAnalysis(config: NetworkConfig):
#     config.outdir = config.outdir /Path(config.analysis_network_mode+'_community_detection')
#     if not config.outdir.exists():
#         config.outdir.mkdir(parents=True, exist_ok=True)
#     # network construction
#     G_raw = network_construction(config.nodefile,config.edgefile)
#     # deal with mode extraction
#     G = network_mode_extraction(G_raw,config.nodefile,config.analysis_network_mode,config.analysis_SEEDNODEID)
#     # Community Detection
#     config.module_info = network_community_detection(G,config.analysis_network_mode,config.comm_detection,config.outdir)
#     # r_script_path = "E:/pro_phos_V2/scripts/3.2functionnal_interaction_network_community_visualisation.R"
#     network_mode = config.visualization_network_mode

#     def to_r_path(path):
#         return str(path).replace("\\", "/")
#     curr_outdir = config.outdir.parent
#     cmd = [
#             "Rscript",
#             to_r_path(str(config.script_path / Path(config.kgbased_functional_network_community_rscript_path))),
#             "--edge_file", to_r_path(curr_outdir/Path('Edge_diff_info.tsv')),
#             "--community_info", to_r_path(curr_outdir / Path('DEP_community_detection') / Path('DEP_community_detection_modules.tsv')),
#             "--outdir", to_r_path(curr_outdir/Path('DEP_community_detection')),
#             "--top_module_num",str(config.kgbased_network_top_n)
#         ]
    
#     if sys.platform.startswith('linux'):
#         conda_base = os.environ.get('CONDA_PREFIX', '/hwfsyt1/SP_MSI/Pipline/software/miniconda3_ppy/miniconda3')
#         # 安全提取 conda root（避免 /envs/ 出现多次）
#         if '/envs/' in conda_base:
#             conda_root = conda_base.split('/envs/')[0]
#         else:
#             conda_root = conda_base  # 可能是 base 环境

#         # 使用 shlex.quote 防止路径含空格或特殊字符
#         import shlex
#         cmd_str = " ".join(shlex.quote(arg) for arg in cmd)
        
#         activate_cmd = (
#             f"source {shlex.quote(conda_root)}/etc/profile.d/conda.sh && "
#             f"conda activate kg_net && "
#             f"{cmd_str}"
#         )
#         full_cmd = ["bash", "-c", activate_cmd]
#     else:
#         # Windows: 直接使用 cmd（假设已通过 conda activate 进入环境）
#         full_cmd = cmd
#         print(">>> DEBUG: 使用了 else 分支（非 Linux）")

#     print(f"📋 执行命令:\n{' '.join(full_cmd)}")
#     result = subprocess.run(full_cmd, capture_output=True, text=True)

#     print("STDOUT:\n", result.stdout)
#     print("STDERR:\n", result.stderr)


# def main():
#     parser = argparse.ArgumentParser(
#         description='对构建所有检测到蛋白和磷酸化蛋白网络进行网络分析',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter
#     )
#     # input files
#     parser.add_argument('-nodefile', required=True,
#                         help='nodes file (TSV)')
#     parser.add_argument('-edgefile', required=True,
#                         help='edges file (TSV)')
#     #the output directory
#     parser.add_argument('-output',type=str,
#                         default=None,
#                         help='指定输出目录（可选），未指定时默认当前目录')
#     #analytical parameter
#     parser.add_argument('-network_mode', type=str,
#                     default='ALL',
#                     choices=["ALL", "KS", "LR", "TF", "DEP", "DESUP","DESDOWN","SEEDNODE"],
#                     help='Select subnetwork type to extract: '
#                          '"ALL" (entire network), '
#                          '"KS" (kinase-substrate interactions), '
#                          '"LR" (ligand-receptor interactions), '
#                          '"TF" (transcription factor-target interactions), '
#                          '"DEP" (differentially expressed genes), '
#                          '"DESUP" (UP regulated genes), '
#                          '"DESDOWN" (DOWN regulated genes), '
#                          '"SEEDNODE" (seed nodes and their neighbors). '
#                          'Default: %(default)s')
#     parser.add_argument('-SEEDNODEID', type=str,
#                     default=None,
#                     help='SEEDNODE ID, Default: %(default)s')
#     parser.add_argument('-comm_detection', type=str,
#                     default='leiden',
#                     choices=["louvain", "leiden", "infomap", "fastgreedy", "walktrap", "LPA", "None"],
#                     help='Community detection algorithm to use: '
#                          '"louvain" (Louvain algorithm), '
#                          '"leiden" (Leiden algorithm), '
#                          '"infomap" (Infomap algorithm), '
#                          '"fastgreedy" (Fast Greedy algorithm), '
#                          '"walktrap" (Walktrap algorithm), '
#                          '"LPA" (Label Propagation Algorithm), '
#                          #######'"spinglass" (Spinglass algorithm). '
#                          'Default: %(default)s (no community detection)')

#     # get all parameters
#     args = parser.parse_args()

#     # parameter validation
#     out_dir = parameter_validation(args.nodefile,args.edgefile,args.output,args.analysis_network_mode)
#     # data loading
#     nodes_df,edges_df = input_load_validation(args.nodefile,args.edgefile,args.network_mode,args.SEEDNODEID)
#     # network construction
#     G_raw = network_construction(nodes_df,edges_df)
#     # deal with mode extraction
#     G = network_mode_extraction(G_raw,nodes_df,args.network_mode,args.SEEDNODEID)
#     # Community Detection
#     network_community_detection(G,args.network_mode,args.comm_detection,out_dir)
    
# if __name__ == "__main__":
#     main()
