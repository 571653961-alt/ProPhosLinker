#!/usr/bin/env python3
"""
ProPhosLinker: Knowledge graph network constructor module (3.2.1).

This module implements the `NetworkConstructor` function (part of ProPhosLinker toolkit), 
which queries Neo4j knowledge graphs to construct protein-phosphoprotein functional 
interaction networks from differential proteomics/phosphoproteomics results.

Core workflow:
1) **Node annotation**: Classifies proteins/phosphosites by:
   - Detection status: proteinDetected/NotproteinDetected
   - Regulation class: UP/DOWN/Non-significant (log2FC ≥ |FC|, FDR < 0.05)
   - Disease association: Queries disease biomarkers from Neo4j
   - Phosphosite mapping: Parent protein → multiple phospho-sites

2) **Edge discovery**: Neo4j Cypher queries for interactions via:
   - `ACTS_ON|CURATED_INTERACTS_WITH` relationships
   - Supports UNIPROTKB or SYMBOL identifiers
   - 14 interaction types (activation, catalysis, phosphorylation, etc.)
   - Deduplicated edge types per source→target pair

3) **Output generation**:
   - `Node_info.tsv`: Comprehensive node metadata (protein + phosphosite info)
   - `Edge_info.tsv`: Source→Target interaction types  
   - `Node_diff_info.tsv` / `Edge_diff_info.tsv`: DE-only subnetwork

Key features:
- **Robust Neo4j integration** with connection pooling and session management
- **Disease-contextualized** network construction 
- **Hierarchical phosphosite representation** (protein → N phospho-sites)
- **Comprehensive validation** of inputs, paths, and parameters
"""

import os
import sys
from pathlib import Path
import pandas as pd
from neo4j import GraphDatabase
from typing import Union, Dict, Literal


class NetworkConfig:
    """
    Configuration class for knowledge graph-based network construction.
    
    Handles input/output paths, omics data parameters, Neo4j connection details,
    and visualization parameters for protein-phosphoprotein networks.
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
        # Identifier type
        identified_type: Literal['UNIPROT', 'SYMBOL'] = 'UNIPROT',
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
        # Network visualization parameters
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
        # Node color scheme
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
            "association": "#d7d6d6",        # Light gray
            "physical association": "#838181",  # Medium gray
            "binding": "#838181",            # Medium gray
            "direct interaction": "#d7d6d6",   # Light gray
            
            "activation": "#d62c0c",         # Red-orange
            "catalysis": "#7f00ff",          # Purple
            "proximity": "#FF9301",          # Orange
            "reaction": "#114335",           # Dark green
            "phosphorylation": "#2FBE95",    # Teal
            "dephosphorylation": "#6B8E23",  # Olive green
            
            "ptmod": "#8C97D6",              # Light purple
            "inhibition": "#0cb6d6",         # Cyan
            "expression": "#FCF402",         # Yellow
            "regulation": "#e4dadd",         # Light purple-gray
            
            "colocalization": "#4c95cd",     # Sky blue
            "covalent binding": "#716F74",   # Violet
            "ubiquitination": "#FF4500",     # Orange-red
            "multiRel": "#9a6728"            # Brown
        },
        # R script paths
        script_path: str = './scripts',
        kgbased_functional_network_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_visualization.R',
        kgbased_functional_network_community_rscript_path: Union[str, Path] = '3.2functionnal_interaction_network_community_visualisation.R',
        function_enrich_rscript_path: Union[str, Path] = '3.2functional_enrichment_function.R',
        # Neo4j connection parameters
        uri: str = "bolt://localhost:7687",
        username: str = "neo4j",
        password: str = "neo4j"
    ):
        # Store all parameters
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
        self.identified_type = identified_type
        self.FC = FC
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
        self.node_up = node_up
        self.node_down = node_down
        self.node_nonsig = node_nonsig
        self.node_notdet = node_notdet
        self.node_disease_border_color = node_disease_border_color
        self.node_notdisease_border_color = node_notdisease_border_color
        self.function_enrich_color_gradient_low = function_enrich_color_gradient_low
        self.function_enrich_color_gradient_high = function_enrich_color_gradient_high
        self.edge_type_colors = edge_type_colors
        self.oudir_default_name = oudir_default_name
        self.script_path = script_path
        self.kgbased_functional_network_rscript_path = kgbased_functional_network_rscript_path
        self.kgbased_functional_network_community_rscript_path = kgbased_functional_network_community_rscript_path
        self.function_enrich_rscript_path = function_enrich_rscript_path
        self.uri = uri
        self.username = username
        self.password = password

        # Validate parameters and load data
        self._validate_params()
        self._data_load()

    def _validate_params(self):
        """Validate input parameters and create output directory."""
        # Convert outdir to Path object and validate
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        elif self.outdir is None:
            self.outdir = Path(os.getcwd())
        elif not isinstance(self.outdir, Path):
            raise TypeError("outdir must be a string or Path object")

        # Check if output directory exists and is writable
        if self.outdir.exists():
            if self.outdir.is_file():
                raise ValueError(f"Output path is a file, not directory: {self.outdir}")
            if not os.access(self.outdir, os.W_OK):
                raise ValueError(f"Output directory exists but is not writable: {self.outdir}")

        # Create default subdirectory
        self.outdir = self.outdir / Path(self.oudir_default_name)
        output_dir = self.outdir
        
        # Ensure output directory exists
        self.outdir.mkdir(parents=True, exist_ok=True)
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f"✅ Output directory ready: {output_dir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {output_dir}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")

        # Validate all input files
        for attr in ['profile', 'phosphoprofile', 'prodiff', 'phosphoprodiff', 'phosphopro_pro_cor']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Path(value))
            elif not isinstance(value, (Path, pd.DataFrame)):
                raise TypeError(f"{attr} must be a string, Path object, or DataFrame")
            
            # Validate file existence and type for Path objects
            if isinstance(value, Path):
                if not value.exists():
                    raise FileNotFoundError(f"{value} does not exist")
                if not value.is_file():
                    raise ValueError(f"{value} is not a valid file")
            
            # Validate DataFrame structure
            if isinstance(value, pd.DataFrame):
                if value.empty:
                    raise ValueError(f"{attr} DataFrame is empty")
                if not all(col in value.columns for col in ['ID', 'Name']):
                    raise ValueError(f"{attr} DataFrame must contain 'ID' and 'Name' columns")

        # Validate fold change threshold
        if self.FC <= 0:
            raise ValueError("Fold change (log2FC) must be > 0")

    def _data_load(self):
        """Load and validate input data files."""
        def datafile_load_validation(file):
            """Load and validate expression profile data."""
            if not isinstance(file, pd.DataFrame):
                # Read tab-separated data with numeric conversion
                df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
                df = df.apply(pd.to_numeric, errors='coerce')
            else:
                df = file
            
            if df.empty:
                raise ValueError(f"{file} expression data is empty")
            if df.index.duplicated().any():
                raise ValueError(f"{file} contains duplicate row indices")
            return df
        
        def difffile_load_validation(file):
            """Load differential expression analysis results."""
            if not isinstance(file, pd.DataFrame):
                df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
            else:
                df = file
            return df

        # Load all input data
        self.profile = datafile_load_validation(self.profile)
        self.phosphoprofile = datafile_load_validation(self.phosphoprofile)
        self.prodiff = difffile_load_validation(self.prodiff)
        self.phosphoprodiff = difffile_load_validation(self.phosphoprodiff)
        
        # Load phosphoprotein-to-protein correspondence
        if not isinstance(self.phosphopro_pro_cor, dict):
            self.phosphopro_pro_cor = pd.read_csv(
                self.phosphopro_pro_cor, sep='\t', header=None, 
                names=['value', 'key']
            ).set_index('key')['value'].to_dict()


def parameter_validation(profile, phosphoprofile, prodiff, phosphoprodiff, 
                        phosphopro_pro_cor, output, FC):
    """
    Legacy parameter validation function for backward compatibility.
    
    Args:
        profile, phosphoprofile, prodiff, phosphoprodiff, phosphopro_pro_cor: Input file paths
        output: Output directory path
        FC: Fold change threshold
        
    Returns:
        Path: Validated output directory path
    """
    # Check input file existence and non-empty status
    def input_file_check(file_path):
        if not os.path.exists(file_path):
            raise FileExistsError(f"File path does not exist: {file_path}")
        if os.path.getsize(file_path) == 0:
            raise ValueError(f"Input file is empty: {file_path}")

    # Validate all input files
    for file_path in [profile, phosphoprofile, prodiff, phosphoprodiff, phosphopro_pro_cor]:
        input_file_check(file_path)

    # Handle output directory
    if output is None:
        base_path = Path(os.getcwd())
    else:
        base_path = Path(output).expanduser().absolute()
        if base_path.exists():
            if base_path.is_file():
                raise ValueError(f"Output path is a file, not directory: {base_path}")
            if not os.access(base_path, os.W_OK):
                raise ValueError(f"Output directory exists but is not writable: {base_path}")

    # Create Network_analysis subdirectory
    default_name = "Network_analysis"
    output_dir = base_path / Path(default_name)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"✅ Output directory ready: {output_dir}", file=sys.stderr)
    except PermissionError:
        raise RuntimeError(f"No permission to create directory: {output_dir}")
    except Exception as e:
        raise RuntimeError(f"Failed to create directory: {str(e)}")

    # Validate fold change
    if FC <= 0:
        raise ValueError("Fold change (log2FC) must be > 0")

    return output_dir


def input_load_validation(profile, phosphoprofile, prodiff, phosphoprodiff, phosphopro_pro_cor):
    """Legacy data loading function for backward compatibility."""
    def datafile_load_validation(file_path):
        df = pd.read_csv(file_path, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
        df = df.apply(pd.to_numeric, errors='coerce')
        if df.empty:
            raise ValueError("Expression data is empty")
        if df.index.duplicated().any():
            raise ValueError("Duplicate row indices found")
        return df
    
    def difffile_load_validation(file_path):
        df = pd.read_csv(file_path, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
        return df

    pro_df = datafile_load_validation(profile)
    phosphopro_df = datafile_load_validation(phosphoprofile)
    prodiff_df = difffile_load_validation(prodiff)
    phosphoprodiff_df = difffile_load_validation(phosphoprodiff)
    phosphopro_pro_cor_dict = pd.read_csv(
        phosphopro_pro_cor, sep='\t', header=None, 
        names=['key', 'value']
    ).set_index('key')['value'].to_dict()
    
    return pro_df, phosphopro_df, prodiff_df, phosphoprodiff_df, phosphopro_pro_cor_dict


def network_constructor(pro_df, phosphopro_df, prodiff_df, phosphoprodiff_df, 
                       phosphopro_pro_cor_dict, FC, identified_type, disease, 
                       uri, username, password):
    """
    Construct protein-phosphoprotein interaction network from knowledge graph.
    
    Args:
        pro_df, phosphopro_df: Protein and phosphoprotein expression data
        prodiff_df, phosphoprodiff_df: Differential expression results
        phosphopro_pro_cor_dict: Phosphosite-to-protein mapping dictionary
        FC: Fold change threshold for significance
        identified_type: 'UNIPROT' or 'SYMBOL' identifier type
        disease: Disease name for association analysis
        uri, username, password: Neo4j connection parameters
        
    Returns:
        tuple: (Node_info dict, Edge_info dict)
    """
    def find_keys_by_value(dictionary, value):
        """Find all keys mapping to a given value in dictionary."""
        return [k for k, v in dictionary.items() if v == value]

    # Extract protein identifiers from data
    protein_names = pro_df.index.tolist()
    phosphoprotein_siteIDs = phosphopro_df.index.tolist()
    phosphoprotein_names = list(set(phosphopro_pro_cor_dict.values()))

    # Classify proteins by differential expression status
    def get_protein_info(diff, FC):
        """Extract fold-change and significance information for proteins."""
        # Extract protein names from site-specific indices
        diff.index = diff.index.to_series().str.split(':').str[0]
        proteinID_FC_Qvalue = {}
        
        for proteinID in diff.index.tolist():
            log2FC = diff.loc[proteinID, 'logFC']
            FDR = diff.loc[proteinID, 'adj.P.Val']
            
            # Classify based on fold-change and FDR thresholds
            if float(log2FC) >= FC and float(FDR) < 0.05:
                CLASS = "UP"
            elif float(log2FC) <= -FC and float(FDR) < 0.05:
                CLASS = "DOWN"
            else:
                CLASS = "Non-significant"
                
            proteinID_FC_Qvalue[proteinID] = {
                'FC': float(log2FC),
                'FDR': float(FDR),
                'CLASS': CLASS
            }
        return proteinID_FC_Qvalue

    proteinID_FC_Qvalue = get_protein_info(prodiff_df, FC)
    phosphoproteinID_FC_Qvalue = get_protein_info(phosphoprodiff_df, FC)

    all_protein_phosphoprotein_names = list(set(protein_names + phosphoprotein_names))

    # Connect to Neo4j database
    driver = GraphDatabase.driver(uri, auth=(username, password))

    # Query disease-associated proteins
    def get_disease_protein_names(identified_type):
        if identified_type == 'UNIPROT':
            cypher_query = """
            MATCH (p:Protein)-[r:IS_BIOMARKER_OF_DISEASE|ASSOCIATED_WITH]->(d:Disease) 
            WHERE toLower(d.name) = toLower($disease)
            RETURN DISTINCT p.id as p_name
            """
        elif identified_type == 'SYMBOL':
            cypher_query = """
            MATCH (p:Protein)-[r:IS_BIOMARKER_OF_DISEASE|ASSOCIATED_WITH]->(d:Disease)
            WHERE toLower(d.name) = toLower($disease)
            RETURN DISTINCT p.name as p_name
            """
        
        disease_proteins = []
        with driver.session() as session:
            result = session.run(cypher_query, disease=disease)
            for record in result:
                disease_proteins.append(record["p_name"])
        return disease_proteins

    disease_protein_names = get_disease_protein_names(identified_type) if disease else []

    # Construct comprehensive node information
    Node_info = {}
    for pro_phospro_name in all_protein_phosphoprotein_names:
        Node_info[pro_phospro_name] = {
            'Pro_type': '',
            'Pro_FC': '',
            'Pro_FDR': '',
            'Pro_class': '',
            'DiseaseRelated': 'yes' if pro_phospro_name in disease_protein_names else 'no',
            'Phosphopro_info': {
                'phosphoSite_detected_num': 0
            }
        }

        # Process detected proteins
        if pro_phospro_name in protein_names and pro_phospro_name in prodiff_df.index:
            Node_info[pro_phospro_name]['Pro_type'] = 'proteinDetected'
            Node_info[pro_phospro_name]['Pro_FC'] = str(proteinID_FC_Qvalue[pro_phospro_name]['FC'])
            Node_info[pro_phospro_name]['Pro_FDR'] = str(proteinID_FC_Qvalue[pro_phospro_name]['FDR'])
            Node_info[pro_phospro_name]['Pro_class'] = proteinID_FC_Qvalue[pro_phospro_name]['CLASS']
            
            if pro_phospro_name in phosphoprotein_names:
                # Multiple phosphosites possible per protein
                phosphoSiteIDs = find_keys_by_value(phosphopro_pro_cor_dict, pro_phospro_name)
                # Filter to significant phosphosites only
                phosphoSiteIDs = [site for site in phosphoSiteIDs 
                                if site in phosphoprodiff_df.index]
                
                if phosphoSiteIDs:
                    Node_info[pro_phospro_name]['Phosphopro_info']['phosphoSite_detected_num'] = len(phosphoSiteIDs)
                    Node_info[pro_phospro_name]['Phosphopro_info']['phosphoSite_proteinIDs'] = phosphoSiteIDs
                    
                    for i, phosphoSiteID in enumerate(phosphoSiteIDs):
                        Node_info[pro_phospro_name]['Phosphopro_info'][phosphoSiteID] = {
                            'phosphoSite_index': i,
                            'Phosphopro_type': 'phosphorylated',
                            'Phosphopro_FC': str(phosphoproteinID_FC_Qvalue[phosphoSiteID]['FC']),
                            'Phosphopro_FDR': str(phosphoproteinID_FC_Qvalue[phosphoSiteID]['FDR']),
                            'Phosphopro_class': phosphoproteinID_FC_Qvalue[phosphoSiteID]['CLASS']
                        }
            else:
                # Protein only, no phosphosites
                Node_info[pro_phospro_name]['Phosphopro_info']['phosphoSite_proteinIDs'] = []
                
        elif pro_phospro_name in phosphoprotein_names:
            # Phosphoprotein only (protein not detected)
            Node_info[pro_phospro_name].update({
                'Pro_type': 'NotproteinDetected',
                'Pro_FC': 'None',
                'Pro_FDR': 'None',
                'Pro_class': 'None'
            })
            
            phosphoSiteIDs = [site for site in find_keys_by_value(phosphopro_pro_cor_dict, pro_phospro_name)
                            if site in phosphoprodiff_df.index]
            
            if phosphoSiteIDs:
                Node_info[pro_phospro_name]['Phosphopro_info'].update({
                    'phosphoSite_detected_num': len(phosphoSiteIDs),
                    'phosphoSite_proteinIDs': phosphoSiteIDs
                })
                
                for i, phosphoSiteID in enumerate(phosphoSiteIDs):
                    Node_info[pro_phospro_name]['Phosphopro_info'][phosphoSiteID] = {
                        'phosphoSite_index': i,
                        'Phosphopro_type': 'phosphorylated',
                        'Phosphopro_FC': str(phosphoproteinID_FC_Qvalue[phosphoSiteID]['FC']),
                        'Phosphopro_FDR': str(phosphoproteinID_FC_Qvalue[phosphoSiteID]['FDR']),
                        'Phosphopro_class': phosphoproteinID_FC_Qvalue[phosphoSiteID]['CLASS']
                    }

    # Extract interaction edges from knowledge graph
    Edge_info = {}
    all_Pro_Nodes = list(Node_info.keys())
    
    for pro_phospro_name in all_Pro_Nodes:
        # Define Cypher query for protein interactions
        if identified_type == 'UNIPROT':
            edge_cypher_query = """
            MATCH (p:Protein {id:$node_id})-[r:ACTS_ON|CURATED_INTERACTS_WITH]-(p2:Protein)
            RETURN DISTINCT
                startNode(r).id as source_name,
                r.interaction_type as rel_type,
                r.action as action,
                endNode(r).id as target_name
            """
            params = {"node_id": pro_phospro_name}
        else:  # SYMBOL
            edge_cypher_query = """
            MATCH (p:Protein {name:$node_id})-[r:ACTS_ON|CURATED_INTERACTS_WITH]-(p2:Protein)
            RETURN DISTINCT
                startNode(r).name as source_name,
                r.interaction_type as rel_type,
                r.action as action,
                endNode(r).name as target_name
            """
            params = {"node_id": pro_phospro_name}
        
        with driver.session() as session:
            result = session.run(edge_cypher_query, **params)
            for record in result:
                source_name = record["source_name"]
                rel_type = record["rel_type"]
                action = record["action"]
                target_name = record["target_name"]
                
                if source_name and target_name and source_name in all_Pro_Nodes and target_name in all_Pro_Nodes:
                    edge_key = f"{source_name}_{target_name}"
                    
                    if edge_key not in Edge_info:
                        Edge_info[edge_key] = []
                    
                    # Add interaction type (prioritize rel_type over action)
                    interaction = rel_type or action
                    if interaction and interaction not in Edge_info[edge_key]:
                        clean_interaction = interaction.replace('phosphorylation reaction', 'phosphorylation')
                        Edge_info[edge_key].append(clean_interaction)

    driver.close()
    return Node_info, Edge_info


def output_save(Node_info, Edge_info, output_dir):
    """
    Save network node and edge information to TSV files.
    
    Args:
        Node_info: Dictionary of node attributes
        Edge_info: Dictionary of edge interactions
        output_dir: Output directory Path object
        
    Returns:
        tuple: (node_df, edge_df) DataFrames for differential network
    """
    # Initialize data lists for DataFrame construction
    f_Node_data = []
    f_Edge_data = []

    # Process node information
    for node_name, node_attrs in Node_info.items():
        pro_type = node_attrs['Pro_type']
        pro_fc = node_attrs['Pro_FC']
        pro_fdr = node_attrs['Pro_FDR']
        pro_class = node_attrs['Pro_class']
        disease_related = node_attrs['DiseaseRelated']
        phospho_info = node_attrs['Phosphopro_info']
        phospho_count = phospho_info['phosphoSite_detected_num']
        
        if phospho_count > 0:
            phospho_sites = phospho_info['phosphoSite_proteinIDs']
            for i, phospho_site in enumerate(phospho_sites):
                phospho_attrs = phospho_info[phospho_site]
                f_Node_data.append([
                    node_name, pro_type, pro_fc, pro_fdr, pro_class, disease_related,
                    phospho_count, i, phospho_site,
                    phospho_attrs['Phosphopro_type'], phospho_attrs['Phosphopro_FC'],
                    phospho_attrs['Phosphopro_FDR'], phospho_attrs['Phosphopro_class']
                ])
        else:
            # No phosphosites detected
            f_Node_data.append([
                node_name, pro_type, pro_fc, pro_fdr, pro_class, disease_related,
                0, 0, '', '', '', '', ''
            ])

    # Create node DataFrame
    node_columns = [
        "NodeName", "Pro_type", "Pro_FC", "Pro_FDR", "Pro_class", "DiseaseRelated",
        "phosphoSite_detected_num", "phosphoSite_index", "phosphoprotein_ID",
        "Phosphopro_type", "Phosphopro_FC", "Phosphopro_FDR", "Phosphopro_class"
    ]
    f_Node_df = pd.DataFrame(f_Node_data, columns=node_columns)

    # Process edge information
    for edge_key, interactions in Edge_info.items():
        # Clean and deduplicate interaction types
        clean_interactions = list(set(';'.join(interactions)
                                    .replace(',', ';')
                                    .replace('binding/association', 'binding')
                                    .split(';')))
        interaction_type = ';'.join(clean_interactions)
        
        if interaction_type:
            source, target = edge_key.split('_')
            f_Edge_data.append([source, target, interaction_type])

    # Create edge DataFrame
    edge_columns = ["Source", "Target", "Interaction_Type"]
    f_Edge_df = pd.DataFrame(f_Edge_data, columns=edge_columns)

    # Save complete network files
    node_path = output_dir / "Node_info.tsv"
    edge_path = output_dir / "Edge_info.tsv"
    f_Node_df.to_csv(node_path, sep='\t', index=False)
    f_Edge_df.to_csv(edge_path, sep='\t', index=False)
    
    print(f"✅ Node and edge information saved: {output_dir}")

    # Extract differential network (UP/DOWN regulated nodes only)
    mask = (
        f_Node_df['Pro_class'].isin(['UP', 'DOWN']) |
        f_Node_df['Phosphopro_class'].isin(['UP', 'DOWN'])
    )
    diff_nodes_tmp = f_Node_df[mask]
    diff_node_names = list(set(diff_nodes_tmp['NodeName'].tolist()))
    
    diff_node_df = f_Node_df[f_Node_df['NodeName'].isin(diff_node_names)]
    diff_edge_df = f_Edge_df[
        (f_Edge_df['Source'].isin(diff_node_names)) & 
        (f_Edge_df['Target'].isin(diff_node_names))
    ]
    
    # Save differential network files
    diff_node_path = output_dir / "Node_diff_info.tsv"
    diff_edge_path = output_dir / "Edge_diff_info.tsv"
    diff_node_df.to_csv(diff_node_path, sep='\t', index=False)
    diff_edge_df.to_csv(diff_edge_path, sep='\t', index=False)
    
    return diff_node_df, diff_edge_df


def NetworkConstructor(config: NetworkConfig):
    """
    Main entry point for network construction using NetworkConfig.
    
    Args:
        config: NetworkConfig instance with all parameters
        
    Returns:
        tuple: (config.nodefile, config.edgefile) populated with DataFrames
    """
    # Construct network from knowledge graph
    Node_info, Edge_info = network_constructor(
        config.profile, config.phosphoprofile, config.prodiff,
        config.phosphoprodiff, config.phosphopro_pro_cor, config.FC,
        config.identified_type, config.disease, config.uri,
        config.username, config.password
    )
    
    # Save results and update config
    config.nodefile, config.edgefile = output_save(Node_info, Edge_info, config.outdir)
    return config.nodefile, config.edgefile
