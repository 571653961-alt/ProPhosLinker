#!/usr/bin/env python3
"""
Differential analysis orchestration module for the ProPhosLinker pipeline.

This module defines the `DifferentialAnalysis` class, which coordinates the full
omics-level differential analysis workflow by wrapping several R-based pipelines
and a Python-based phosphorylation-rate subtyping step.

Core responsibilities:
- Build and execute the omics differential analysis (step 2.1) via an R script:
  * Takes protein and phosphoprotein matrices, mapping file, and group information.
  * Applies log2 fold-change and p-value thresholds.
  * Produces differential result tables and standard visualizations.
- Run phosphorylation-rate subtyping (step 2.2) in Python:
  * Loads differential results from step 2.1.
  * Computes phosphorylation-rate–based patterns and subtypes using quantile-based
    logic and clustering.
- Run the stable differential network analysis (step 2.3) via an R script:
  * Uses differential results to build a stable differential network.
  * Applies multiple stability and edge-filtering parameters.
  * Generates network-level outputs and plots.

Environment handling:
- On Linux, the omics differential analysis can be executed inside a dedicated
  conda environment (e.g. `differential_analysis`) if configured.
- On other platforms, the module calls `Rscript` directly and assumes it is
  available on the system PATH.

Typical usage:
- Instantiate `DifferentialAnalysis` with a fully populated `DifferentialConfig`
  object, then call `run()` to execute steps 2.1 → 2.2 → 2.3 as a single pipeline.
"""


from ..config import DifferentialConfig, ResultDirConfig, ScriptConfig
import subprocess
from pathlib import Path


def to_r_path(path: str) -> str:
    """Convert a local path to an R-friendly path string."""
    return str(Path(path)).replace('\\', '/')

class DifferentialAnalysis:
    """
    Differential analysis pipeline runner.

    Notes
    - This class is a thin orchestrator around R-based pipelines.
    - On Linux, the implementation optionally activates specific conda
      environments before calling `Rscript`.
    - On non-Linux platforms, it calls `Rscript` directly (expects it on PATH).
    """
    def __init__(self, config: DifferentialConfig):
        self.config = config
    
    def omics_differential_analysis(self):
        """
        2.1 Omics differential analysis.

        Runs the R workflow to compute differential signals for protein and
        phosphosite matrices and generate standard visualizations/summary tables.
        """
        print("📊 Starting omics differential analysis...")
        
        outdir = Path(self.config.outdir).parent / Path(ResultDirConfig.omics_differential_analysis)    #Path("2.1Omics_Differential_Analysis")
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "Rscript", to_r_path(ScriptConfig.get_full_path('omics_differential_rscript_path')),
            "--profile", to_r_path(self.config.profile),
            "--phosfile", to_r_path(self.config.phosfile),
            "--pro_phos_cor", to_r_path(self.config.mappingfile),
            "--compare_groups", to_r_path(self.config.sample_group),
            "--outdir", to_r_path(outdir),
            "--group_comparing", self.config.group_comparing,
            "--log2FC_pro", str(self.config.pro_log2FC),
            "--log2FC_phos", str(self.config.phos_log2FC),
            "--p_val", str(self.config.diff_p_val),
            "--omics1_name", self.config.omics1_name,
            "--omics2_name", self.config.omics2_name,
            "--quadrant_plot_up_color", self.config.omics_diff_quadrant_pro_color,
            "--quadrant_plot_down_color", self.config.omics_diff_quadrant_phos_color,
            # "--verbose"
        ]

        # Conditionally add --prodiff
        if self.config.prodiff:
            cmd.extend(["--prodiff", to_r_path(self.config.prodiff)])

        # Conditionally add --phosdiff
        if self.config.phosdiff:
            cmd.extend(["--phosdiff", to_r_path(self.config.phosdiff)])


        
        full_cmd = cmd

        print(f"📋 Executing command:\n{' '.join(full_cmd)}")
        result = subprocess.run(full_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Omics differential analysis completed")
            if result.stdout:
                print(f"📄 Output:\n{result.stdout}")
            return True
        else:
            print("❌ Omics differential analysis failed")
            print(f"📄 STDOUT:\n{result.stdout}")
            print(f"❌ STDERR:\n{result.stderr}")
            return False
        return True
    
    def phosrate_subtyping(self):
        """
        2.2 Phosphorylation-rate subtyping.

        Builds a `PhosRateSubtyping` configuration and runs the quantile-based
        phosphorylation-rate analysis pipeline.
        """
        print("🧬 Starting phosphorylation-rate subtyping...")
        
        try:
            # Lazy import to keep the base import surface light
            from . import phos_rate_subtyping
            
            # Expected filenames produced by step 2.1
            group1, group2 = self.config.group_comparing.split(':')
            prodiff_file = f"{group1}_vs_{group2}_{self.config.omics1_name}.tsv"
            phosphoprodiff_file = f"{group1}_vs_{group2}_{self.config.omics2_name}.tsv"
            
            # Build configuration object
            PhosratesubtypingConfig = phos_rate_subtyping.PhosRateSubtyping(
                profile=self.config.profile,
                phosphoprofile=self.config.phosfile,
                phosphopro_pro_cor=self.config.mappingfile,
                prodiff= Path(self.config.outdir).parent / Path(ResultDirConfig.omics_differential_analysis) / Path(prodiff_file),
                phosphoprodiff= Path(self.config.outdir).parent / Path(ResultDirConfig.omics_differential_analysis) / Path(phosphoprodiff_file),
                group_info=self.config.sample_group,
                outdir=Path(self.config.outdir).parent / Path(ResultDirConfig.phosrate_subtyping),
                omics1_name=self.config.omics1_name,
                omics2_name=self.config.omics2_name,
                group_vs=self.config.group_comparing,
                phos_rate_FC=self.config.phos_rate_FC,
                quantile_rank=self.config.quantile_rank,
                quantile_rank_step=self.config.quantile_rank_step,
                plot_sites_num=self.config.plot_sites_num,
                frac_param=self.config.frac_param,
                cluster_num=self.config.cluster_num,
                vs_total_phos_rate_upper=self.config.vs_total_phos_rate_upper,
                vs_total_phos_rate_lower=self.config.vs_total_phos_rate_lower,
                top_tail_diff_site_num=self.config.top_tail_diff_site_num,
                r_script_path=ScriptConfig.get_full_path('phosRate_quantile_subtyping_rscript_path'),
                raw_pro_point_color=self.config.raw_pro_point_color,
                fitted_pro_point_color=self.config.fitted_pro_point_color,
                fitted_pro_line_color=self.config.fitted_pro_line_color,
                pro_fc_line_color=self.config.pro_fc_line_color,
                total_phos_rate_line_color=self.config.total_phos_rate_line_color,
                color_lst=self.config.color_lst,
                Mfuzz_pluscolor=self.config.Mfuzz_pluscolor,
                Mfuzz_minuscolor=self.config.Mfuzz_minuscolor,
                Mfuzz_allcolor=self.config.Mfuzz_allcolor,
                Mfuzz_onlyupcolor=self.config.Mfuzz_onlyupcolor,
                Mfuzz_onlydowncolor=self.config.Mfuzz_onlydowncolor,
                Mfuzz_updowncolor=self.config.Mfuzz_updowncolor,
                Mfuzz_notsigcolor=self.config.Mfuzz_notsigcolor
            )
            
            # 执行磷酸化率分型分析
            success = phos_rate_subtyping.phos_rate_subtyping(PhosratesubtypingConfig)
          
            if success:
                print("✅ Phosphorylation-rate subtyping completed")
            else:
                print("❌ Phosphorylation-rate subtyping failed")
            
            return success
            
        except ImportError:
            print("❌ Failed to import `phos_rate_subtyping` module")
            return False
        except Exception as e:
            print(f"❌ Phosphorylation-rate subtyping error: {e}")
            return False
    
    def stable_differential_network(self):
        """
        2.3 Stable differential network analysis.

        Runs the R pipeline to construct a stable differential network and
        produce network-level outputs/visualizations.
        """
        print("🔗 Starting stable differential network analysis...")
        
        outdir = Path(self.config.outdir).parent / Path(ResultDirConfig.stable_differential_network)
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)

        # Build command: pass raw values (do not manually quote tokens)
        cmd = [
            "Rscript", to_r_path(ScriptConfig.get_full_path('differential_network_rscript_path')),
            "--script_path", to_r_path(Path(self.config.script_path) / '2.3Stable_Differential_Network'),
            "--profile", to_r_path(self.config.profile),
            "--phosfile", to_r_path(self.config.phosfile),
            "--pro_phos_cor", to_r_path(self.config.mappingfile),
            "--sample_group", to_r_path(self.config.sample_group),
            "--diff_pro_path", to_r_path(self.config.differential_network_diff_pro_path),
            "--diff_phos_path", to_r_path(self.config.differential_network_diff_phos_path),
            "--outdir", to_r_path(outdir),
            "--group_comparing", self.config.group_comparing,
            "--omics1_name", self.config.omics1_name,
            "--omics2_name", self.config.omics2_name,
            "--filter_num", str(self.config.differential_network_filter_num),
            "--FC_threshold", str(self.config.differential_network_FC_threshold),
            "--p_threshold", str(self.config.differential_network_p_threshold),
            "--p_value_type", self.config.differential_network_p_value_type,
            "--nBoots", str(self.config.differential_network_nBoots),
            "--bootnet_R_threshold", str(self.config.differential_network_bootnet_R_threshold),
            "--nCores", str(self.config.differential_network_nCores),
            "--stability_threshold", str(self.config.differential_network_stability_threshold),
            "--cor_method", self.config.differential_network_cor_method,
            "--edge_FC_threshold", str(self.config.differential_network_edge_FC_threshold),
            "--edge_p_threshold", str(self.config.differential_network_edge_p_threshold),
            "--max_subnet_num", str(self.config.differential_network_max_subnet_num),
            "--enrich_fromType", str(self.config.identified_type),
            "--R_threshold", str(self.config.differential_network_R_threshold),
            "--edge_color_pos", self.config.differential_network_edge_color_pos,
            "--edge_color_neg", self.config.differential_network_edge_color_neg,
            "--Enhanced_in_N", self.config.differential_network_Enhanced_in_N,
            "--Enhanced_in_T", self.config.differential_network_Enhanced_in_T,
            "--Only_in_N", self.config.differential_network_Only_in_N,
            "--Only_in_T", self.config.differential_network_Only_in_T,
            "--Conflict_relation", self.config.differential_network_Conflict_relation,
            "--fill_gradientn_color", ";".join(self.config.differential_network_gradientn_color),
            "--color_gradient_low", self.config.differential_network_function_enrichment_color_gradient_low,
            "--color_gradient_high", self.config.differential_network_function_enrichment_color_gradient_high,
            # "--verbose"
        ]
        
        full_cmd = cmd
        print(f"📋 Executing command:\n{' '.join(str(arg) for arg in full_cmd)}")
        result = subprocess.run(full_cmd, capture_output=True, text=True)

        # Handle results
        
        if result.returncode == 0:
            print("✅ Stable differential network analysis completed")
            if result.stdout:
                print(f"📄 Output:\n{result.stdout}")
            return True
        else:
            print("❌ Stable differential network analysis failed")
            print(f"📄 STDOUT:\n{result.stdout}")
            print(f"❌ STDERR:\n{result.stderr}")
            return False
    
    def validate_config(self):
        """Validate required configuration fields and input file existence."""
        required_attributes = [
            'script_path', 'profile', 'phosfile', 'mappingfile', 'sample_group',
            'outdir', 'omics1_name', 'omics2_name', 'group_comparing'
        ]
        
        for attr in required_attributes:
            if not hasattr(self.config, attr) or getattr(self.config, attr) is None:
                raise ValueError(f"Missing required configuration: {attr}")
        
        # Validate input files (paths) exist when they are provided as strings/Paths
        for file_attr in ['profile', 'phosfile', 'mappingfile', 'sample_group']:
            file_path = getattr(self.config, file_attr)
            if isinstance(file_path, (str, Path)) and not Path(file_path).exists():
                raise FileNotFoundError(f"File not found for {file_attr}: {file_path}")
    
    def run(self):
        """Run the full differential analysis pipeline (2.1 → 2.2 → 2.3)."""
        try:
            print("📈 Starting differential analysis pipeline...")
            
            # 1) Validate config
            self.validate_config()
            
            # 2) Run sub-analyses
            success_omics = self.omics_differential_analysis()
            success_phosrate = self.phosrate_subtyping()
            success_network = self.stable_differential_network()
            
            # 3) Summarize
            overall_success = success_omics and success_phosrate and success_network
            
            if overall_success:
                print("🎉 Differential analysis pipeline completed")
            else:
                print("⚠️  Differential analysis pipeline partially completed")
                if not success_omics:
                    print("  - Omics differential analysis failed")
                if not success_phosrate:
                    print("  - Phosphorylation-rate subtyping failed")
                if not success_network:
                    print("  - Stable differential network analysis failed")
            
            return overall_success
            
        except Exception as e:
            print(f"❌ Differential analysis pipeline error: {e}")
            return False