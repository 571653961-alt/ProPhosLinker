#!/usr/bin/env python3
"""
ProPhosLinker: Functional analysis orchestration module for the pro_phos_cor pipeline.

This module defines the `FunctionalAnalysis` class (part of ProPhosLinker toolkit), 
which coordinates the complete functional interpretation workflow using differential
analysis results from previous pipeline stages.

Core responsibilities (executes steps 3.1 → 3.2 → 3.3):
1) **Omics functional enrichment** (3.1): R-based GO/KEGG enrichment analysis on
   protein and phosphoprotein differential results with customizable p-value and
   category display thresholds.

2) **Knowledge graph functional network** (3.2): Builds disease-contextualized 
   functional interaction networks using Neo4j-backed knowledge graph data, with
   configurable community detection, node filtering, and visualization parameters.

3) **Hub protein expansion network** (3.3): Expands core hub proteins from the
   knowledge graph network into comprehensive functional networks with enrichment
   visualization and customizable layout/topology parameters.

Key features:
- Direct `Rscript` execution for enrichment analysis (assumes R on PATH).
- Lazy imports for Python-based network modules to minimize import surface.
- Comprehensive input validation and dependency checking.
- Cross-platform path handling via `to_r_path()` utility.
- Detailed success/failure reporting for each sub-analysis.

Typical usage:
- Construct `FunctionalConfig` with required paths and Neo4j credentials.
- Create `FunctionalAnalysis(config)` and call `run()` to execute the full
  functional analysis workflow.
"""


from ..config import FunctionalConfig,ResultDirConfig,ScriptConfig
import subprocess
from pathlib import Path


def to_r_path(path: str) -> str:
    """
    Convert a filesystem path to an R-friendly format by replacing backslashes
    with forward slashes.
    """
    return str(Path(path)).replace('\\', '/')


class FunctionalAnalysis:
    def __init__(self, config: FunctionalConfig):
        self.config = config
    
    def omics_function_enrichment(self):
        """3.1 Omics functional enrichment analysis.

        Runs an R script to perform GO/KEGG enrichment analysis on proteomics and
        phosphoproteomics differential expression results.
        """
        # Status message: starting omics functional enrichment analysis
        print("🧬 Starting omics functional enrichment analysis...")
        
        # Prepare the output directory for omics functional enrichment results.
        # The results are stored in a subdirectory of the main output directory.
        outdir = Path(self.config.outdir).parent / Path(ResultDirConfig.omics_function_enrichment)
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)
        
        # Build the command line arguments that will be passed to the R script.
        # We convert all paths to R-friendly format and forward configuration
        # parameters from the FunctionalConfig.
        cmd = [
            "Rscript", to_r_path(ScriptConfig.get_full_path('omics_function_enrichment_rscript_path')),
            "--pro_diff", to_r_path(self.config.prodiff),
            "--phos_diff", to_r_path(self.config.phosdiff),
            "--phos_pro", to_r_path(self.config.mappingfile),
            "--outdir", to_r_path(outdir),
            "--log2FC", str(self.config.network_log2FC),
            "--diff_p_adj", str(self.config.diff_p_adj),
            "--GO_showCategory", str(self.config.GO_showCategory),
            "--KEGG_showCategory", str(self.config.KEGG_showCategory),
            "--pvalueCutoff", str(self.config.pvalueCutoff),
            "--omics1_name", self.config.omics1_name,
            "--omics2_name", self.config.omics2_name,
            "--enrich_fromType", self.config.identified_type,
            "--color_gradient_low", self.config.omics_function_enrichment_color_gradient_low,
            "--color_gradient_high", self.config.omics_function_enrichment_color_gradient_high,
            # "--verbose"
        ]


        full_cmd = cmd

        print(f"📋 Running command:\n{' '.join(full_cmd)}")
        result = subprocess.run(full_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Omics functional enrichment analysis completed")
            if result.stdout:
                print(f"📄 Output:\n{result.stdout}")
            return True
        else:
            print("❌ Omics functional enrichment analysis failed")
            print(f"📄 Stdout:\n{result.stdout}")
            print(f"❌ Stderr:\n{result.stderr}")
            return False
    
    def kgbased_functional_network(self):
        """3.2 Knowledge graph-based functional network analysis.

        This step builds a functional interaction network based on knowledge graph data
        and the differential analysis results, then visualizes the network.
        """
        # Status message: starting knowledge graph-based functional network analysis
        print("🔗 Starting knowledge graph functional network analysis...")
        
        try:
            # Dynamically import the functionnal_interaction_network module
            from . import functionnal_interaction_network
            
            # Build a configuration object for the functional interaction network analysis.
            # This encapsulates all parameters required by the downstream module.
            network_config = functionnal_interaction_network.NetworkConfig(
                script_path=self.config.script_path,
                profile=self.config.profile,
                phosphoprofile=self.config.phosfile,
                nodefile=None,
                edgefile=None,
                prodiff=self.config.prodiff,
                phosphoprodiff=self.config.phosdiff,
                phosphopro_pro_cor=self.config.mappingfile,
                outdir= Path(self.config.outdir).parent / Path(ResultDirConfig.kgbased_functional_network),
                omics1_name=self.config.omics1_name,
                omics2_name=self.config.omics2_name,
                identified_type=self.config.identified_type,
                kgbased_network_top_n=self.config.kgbased_network_top_n,
                FC=self.config.network_log2FC,
                disease=self.config.kgbased_network_disease,
                analysis_network_mode='DEP',
                analysis_SEEDNODEID=None,
                comm_detection = self.config.kgbased_network_comm_detection,
                top_nodes_visualization_num = self.config.kgbased_network_top_nodes_visualization_num,
                module_id=0,
                max_phosphoSite_displayed = self.config.kgbased_network_max_phosphoSite_displayed,
                visualization_network_mode='ALL',
                visualization_SEEDNODEID=None,
                node_filtering = self.config.kgbased_network_node_filtering,
                network_layout = self.config.kgbased_network_layout,
                #
                uri=self.config.uri,
                username=self.config.username,
                password=self.config.password,
                node_up=self.config.kgbased_functional_network_node_up,
                node_down=self.config.kgbased_functional_network_node_down,
                node_nonsig=self.config.kgbased_functional_network_node_nonsig,
                node_notdet=self.config.kgbased_functional_network_node_notdet,
                node_disease_border_color=self.config.kgbased_functional_network_node_disease_border_color,
                node_notdisease_border_color=self.config.kgbased_functional_network_node_notdisease_border_color,
                edge_type_colors=self.config.kgbased_functional_network_edge_type_colors,
                function_enrich_color_gradient_low=self.config.kgbased_functional_network_function_enrich_color_gradient_low,
                function_enrich_color_gradient_high=self.config.kgbased_functional_network_function_enrich_color_gradient_high
            )
            # Execute the functional network analysis
            success = functionnal_interaction_network.functionnal_interaction_network(network_config)
            
            if success:
                print("✅ Knowledge graph functional network analysis completed")
            else:
                print("❌ Knowledge graph functional network analysis failed")
            
            return success
            
        except ImportError:
            print("❌ Could not import functionnal_interaction_network module")
            return False
        except Exception as e:
            print(f"❌ Error running knowledge graph functional network analysis: {e}")
            return False
    
    def hub_protein_functional_network(self):
        """3.3 Hub protein functional network analysis.

        This step builds and visualizes an expanded hub protein network based on
        the knowledge graph functional network results.
        """
        # Status message: starting hub protein functional network analysis
        print("⭐ Starting hub protein functional network analysis...")
        
        try:
            # Dynamically import the hub_protein_visualization module
            from . import hub_protein_visualization
            
            # Construct expected input paths for node/edge differential information.
            # These files are produced by the knowledge graph functional network step.
            node_path = Path(self.config.outdir).parent / Path(ResultDirConfig.kgbased_functional_network) / Path('Node_diff_info.tsv')
            edge_path = Path(self.config.outdir).parent / Path(ResultDirConfig.kgbased_functional_network) / Path('Edge_diff_info.tsv')
            output_dir = Path(self.config.outdir).parent / Path(ResultDirConfig.hub_protein_functional_network)
            
            # Run the hub protein expansion and visualization module.
            success = hub_protein_visualization.hub_protein_expansion(
                node_path=node_path,
                edge_path=edge_path,
                outdir_parent=output_dir,
                top_n=self.config.hub_protein_top_n,
                d_num=self.config.hub_protein_d_num,
                network_mode=self.config.hub_protein_network_mode,
                network_layout=self.config.hub_protein_network_layout,
                r_script_path=ScriptConfig.get_full_path('hub_protein_functional_network_rscript_path'),
                omics1_name=self.config.omics1_name,
                omics2_name=self.config.omics2_name,
                identified_type = self.config.identified_type,
                function_enrich_rscript_path=ScriptConfig.get_full_path('function_enrich_rscript_path'),
                hub_protein_network_node_up=self.config.hub_protein_network_node_up,
                hub_protein_network_node_down=self.config.hub_protein_network_node_down,
                hub_protein_network_node_nonsig=self.config.hub_protein_network_node_nonsig,
                hub_protein_network_node_notdet=self.config.hub_protein_network_node_notdet,
                hub_protein_network_node_disease_border_color=self.config.hub_protein_network_node_disease_border_color,
                hub_protein_network_node_notdisease_border_color=self.config.hub_protein_network_node_notdisease_border_color,
                hub_protein_network_function_enrich_color_gradient_low=self.config.hub_protein_network_function_enrich_color_gradient_low,
                hub_protein_network_function_enrich_color_gradient_high=self.config.hub_protein_network_function_enrich_color_gradient_high,
                hub_protein_network_edge_type_colors=self.config.hub_protein_network_edge_type_colors
            )
                
            if success:
                print("✅ Hub protein functional network analysis completed")
            else:
                print("❌ Hub protein functional network analysis failed")
            
            return success
            
        except ImportError:
            print("❌ Could not import hub_protein_visualization module")
            return False
        except Exception as e:
            print(f"❌ Error running hub protein functional network analysis: {e}")
            return False
    
    def validate_config(self):
        """Validate required configuration attributes and input files.

        Ensures that all required configuration fields are present and that
        the expected input files exist on disk.
        """
        required_attributes = [
            'script_path', 'profile', 'phosfile', 'mappingfile', 
            'prodiff', 'phosdiff', 'outdir', 'omics1_name', 'omics2_name'
        ]
        
        for attr in required_attributes:
            if not hasattr(self.config, attr) or getattr(self.config, attr) is None:
                raise ValueError(f"Missing required configuration: {attr}")
        
        # Verify input files exist
        for file_attr in ['profile', 'phosfile', 'mappingfile', 'prodiff', 'phosdiff']:
            file_path = getattr(self.config, file_attr)
            if isinstance(file_path, (str, Path)) and not Path(file_path).exists():
                raise FileNotFoundError(f"File not found for {file_attr}: {file_path}")
    
    def check_dependencies(self):
        """Check that required internal Python modules can be imported.

        This helps detect missing dependencies before running the full analysis.
        """
        try:
            from . import functionnal_interaction_network
            from . import hub_protein_visualization
            # Status message: all required internal modules are installed
            print("✅ All required internal modules are available")
            return True
        except ImportError as e:
            # Status message: missing a required module
            print(f"❌ Missing dependency module: {e}")
            return False
    
    def run(self):
        """Run the full functional analysis workflow.

        Validates configuration, checks dependencies, and executes the three
        main analysis stages in order.
        """
        try:
            # Status message: starting the full functional analysis workflow
            print("🧬 Starting functional analysis workflow...")
            
            # 1. Validate configuration
            self.validate_config()
            
            # 2. Check dependencies
            if not self.check_dependencies():
                return False
            
            # 3. Run the three sub-analyses
            success_enrichment = self.omics_function_enrichment()
            success_kg_network = self.kgbased_functional_network()
            success_hub_network = self.hub_protein_functional_network()
            
            # 4. Aggregate results
            overall_success = success_enrichment and success_kg_network and success_hub_network
            
            if overall_success:
                print("🎉 Functional analysis workflow completed")
                print(f"📁 Results saved to: {self.config.outdir}")
            else:
                print("⚠️ Functional analysis workflow partially completed")
                if not success_enrichment:
                    print("  - Omics functional enrichment failed")
                if not success_kg_network:
                    print("  - Knowledge graph functional network analysis failed")
                if not success_hub_network:
                    print("  - Hub protein functional network analysis failed")
            
            return overall_success
            
        except Exception as e:
            print(f"❌ Error running functional analysis workflow: {e}")
            return False