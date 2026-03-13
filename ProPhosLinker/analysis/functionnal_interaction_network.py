#!/usr/bin/env python3
"""
ProPhosLinker: Knowledge graph-based functional interaction network analysis module (step 3.2).

This module defines the `NetworkConfig` class and `functionnal_interaction_network()` orchestrator 
(part of ProPhosLinker toolkit), providing comprehensive protein-phosphoprotein interaction 
network construction, analysis, and visualization using Neo4j-backed knowledge graphs.

Core workflow (3-stage pipeline):
1) **Network Construction** (`pro_phosphopro_network_contructor.NetworkConstructor`):
   - Filters DEPs/DEP sites by log2FC threshold and disease context
   - Maps phosphosites to parent proteins via correspondence table
   - Queries Neo4j for functional interactions (14 edge types: activation, catalysis, phosphorylation, etc.)

2) **Network Analysis** (`pro_phosphopro_network_analysis.NetworkAnalysis`):
   - **Modes**: ALL, KS, LR, TF, DEP(UP/DOWN), SEEDNODE
   - **Community detection**: louvain, leiden, infomap, fastgreedy, walktrap, LPA
   - **Node centrality**: betweenness, degree, closeness, eigenvector, pagerank, etc.

3) **Network Visualization** (`pro_phosphopro_network_visualization.NetworkVisualization`):
   - **Layouts**: force-directed (fr/kk/dh), tree, circular, grid
   - **Node styling**: 4 regulation states (UP/DOWN/NONSIG/NOTDET) + disease borders
   - **Edge coloring**: 14 interaction types with semantic color scheme
   - **Controls**: top-N nodes, max phospho-sites per protein, modular views

Advanced features:
- **Dual ID support**: UNIPROTKB or SYMBOL mapping
- **Rich Neo4j integration**: bolt:// protocol with configurable credentials
- **Modular R visualization**: Separate scripts for main network and community views
- **Comprehensive validation**: Input files, paths, parameters, write permissions

Typical usage:
- Called from `FunctionalAnalysis.kgbased_functional_network()` after omics enrichment (step 3.1)
- Generates `Node_diff_info.tsv`, `Edge_diff_info.tsv` for downstream hub expansion (step 3.3)
"""


import pandas as pd
from typing import Union, Dict,Literal
from pathlib import Path
import os
import sys


class NetworkConfig:
    def __init__(
        self,
        *,
        # Input files
        profile : Union[str, Path, pd.DataFrame],
        phosphoprofile : Union[str, Path, pd.DataFrame],
        prodiff : Union[str, Path, pd.DataFrame],
        phosphoprodiff : Union[str, Path, pd.DataFrame],
        phosphopro_pro_cor : Union[str, Path, Dict],
        # Output directory
        outdir : Union[str, Path],
        # Parameters related to network data
        nodefile : Union[str, Path, pd.DataFrame] = None,
        edgefile : Union[str, Path, pd.DataFrame] = None,
        # 
        module_info : Union[str, Path, pd.DataFrame] = None,
        #omics name
        omics1_name : str = 'Pro',
        omics2_name : str = 'Phos',
        #
        identified_type : Literal['UNIPROT', 'SYMBOL'] = 'UNIPROT',
        # the parameter of the network constructor 
        FC : float = 1.5,
        disease : str = '',
        kgbased_network_top_n : int = 6,
        # the parameter of the network analysis
        analysis_network_mode : Literal["ALL",    # 3 catalogs
                        "KS","LR","TF",
                        "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'DEP',
        analysis_SEEDNODEID : str = None,
        comm_detection : Literal["louvain","leiden","infomap","fastgreedy","walktrap", "LPA"] ='infomap',
        #  the parameter of the network visualization
        top_nodes_visualization_num: int = 20,
        module_id: int = 0,
        max_phosphoSite_displayed: int = 5,
        visualization_network_mode : Literal["ALL", 
                        "KS","LR","TF",
                        "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'ALL',
        visualization_SEEDNODEID : str = None,
        node_filtering : Literal["betweenness", "degree",  "closeness", "harmonic","eigenvector", "pagerank", "alpha", "hub",  "authority"] = 'betweenness',
        network_layout : Literal["fr","kk","dh","stress","tree","gem","graphopt", "lgl","circle","grid"] = 'kk',
        node_up :str = "#B83A2D",
        node_down :str = "#657f68",
        node_nonsig :str = "#E3C79F",
        node_notdet :str = "grey90",
        node_disease_border_color :str = "#535c54",
        node_notdisease_border_color :str = "#acabab",
        function_enrich_color_gradient_low :str = "#175663",
        function_enrich_color_gradient_high :str = "#90362d",
        # Edge type color mapping (for graph visualization)
        edge_type_colors: Dict[str, str] = {
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
        },
        #
        script_path : str = './scripts',
        kgbased_functional_network_rscript_path : Union[str, Path] = '3.2functionnal_interaction_network_visualization.R',
        kgbased_functional_network_community_rscript_path : Union[str, Path] = '3.2functionnal_interaction_network_community_visualisation.R',
        function_enrich_rscript_path : Union[str, Path] = '3.2functional_enrichment_function.R',
        
        ####
        uri : str =  "bolt://localhost:7687",
        username : str =  "neo4j",
        password : str =  "neo4j"
        ):
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
        #
        self.identified_type = identified_type
        #
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
        #
        self.node_up = node_up
        self.node_down = node_down
        self.node_nonsig = node_nonsig
        self.node_notdet = node_notdet
        self.node_disease_border_color = node_disease_border_color
        self.node_notdisease_border_color = node_notdisease_border_color
        self.function_enrich_color_gradient_low = function_enrich_color_gradient_low
        self.function_enrich_color_gradient_high = function_enrich_color_gradient_high
        # Edge type color mapping
        self.edge_type_colors = edge_type_colors
        #
        self.script_path = script_path
        self.kgbased_functional_network_rscript_path  = kgbased_functional_network_rscript_path
        self.kgbased_functional_network_community_rscript_path  = kgbased_functional_network_community_rscript_path
        self.function_enrich_rscript_path = function_enrich_rscript_path
        #
        self.uri = uri
        self.username = username
        self.password = password

        self._validate_params()
        self._data_load()

    def _validate_params(self):
        # Ensure outdir is a Path object
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        elif self.outdir is None:
            self.outdir = os.getcwd()
        elif not isinstance(self.outdir, Path):
            raise TypeError("outdir must be a string or Path object")
        if self.outdir.exists():
            if self.outdir.is_file():
                raise ValueError(f"Output path is a file, not a directory: {self.outdir}")
            if not os.access(self.outdir, os.W_OK):
                raise ValueError(f"Output directory exists and is not writable: {self.outdir}")
        # Create the output directory if it does not exist
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        try:
            self.outdir.mkdir(parents=True,exist_ok=True)
            # Status message: output directory is ready
            print(f"✅ Output directory ready: {self.outdir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {self.outdir}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")
        
        # Ensure all input files are Path objects or DataFrames
        for attr in ['profile', 'phosphoprofile', 'prodiff', 'phosphoprodiff', 'phosphopro_pro_cor']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Path(value))
            elif not isinstance(value, (Path, pd.DataFrame)):
                raise TypeError(f"{attr} must be a string, Path object, or DataFrame")
            # Ensure the input files exist if they are Path objects
            if isinstance(value, Path) and not value.exists():
                raise FileNotFoundError(f"{value} does not exist")
            # Ensure the input files are readable if they are Path objects
            if isinstance(value, Path) and not value.is_file():
                raise ValueError(f"{value} is not a valid file")
            # Ensure the input files are DataFrames if they are not Path objects
            if isinstance(value, pd.DataFrame):
                if value.empty:
                    raise ValueError(f"{attr} DataFrame is empty")
                if not all(col in value.columns for col in ['ID', 'Name']):
                    raise ValueError(f"{attr} DataFrame must contain 'ID' and 'Name' columns")
    
        if self.FC <= 0:
            raise ValueError("FC (log2FC) must be > 0")

    def _data_load(self):
        def datafile_load_validation(file):
            if not isinstance(file, pd.DataFrame):
                    # Load data from tab-delimited file
                    df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
                    df = df.apply(pd.to_numeric, errors='coerce')
            else:
                df = file
            if df.empty:
                raise ValueError(file.as_posix() + " expression profile data is empty")
            if df.index.duplicated().any():
                raise ValueError(file.as_posix() + " has duplicated row indices")
            return df
        
        def difffile_load_validation(file):
            if not isinstance(file, pd.DataFrame):
                    df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
            else:
                df = file
            return df
        
        self.profile = datafile_load_validation(self.profile)
        self.phosphoprofile = datafile_load_validation(self.phosphoprofile)
        self.prodiff = difffile_load_validation(self.prodiff)
        self.phosphoprodiff = difffile_load_validation(self.phosphoprodiff)
        if not isinstance(self.phosphopro_pro_cor, dict):
            self.phosphopro_pro_cor = pd.read_csv(self.phosphopro_pro_cor,sep='\t',header=None,names=['value','key']).set_index('key')['value'].to_dict()




def functionnal_interaction_network(config):
    try:
        from . import pro_phosphopro_network_contructor,pro_phosphopro_network_analysis,pro_phosphopro_network_visualization
        pro_phosphopro_network_contructor.NetworkConstructor(config)
        pro_phosphopro_network_analysis.NetworkAnalysis(config)
        pro_phosphopro_network_visualization.NetworkVisualization(config)
        return True
    except Exception as e:
        print(f"functionnal_interaction_network error: {e}", file=sys.stderr)
        return False

