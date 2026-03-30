#!/usr/bin/env python3
"""
ProPhosLinker: Pattern analysis orchestration module for the ProPhosLinker pipeline.

This module defines the `PatternAnalysis` class (part of ProPhosLinker toolkit), 
which wraps an R-based pattern analysis workflow driven by a `PatternAnalysisConfig` object.

Core responsibilities:
- Build an `Rscript` command from the configuration and call the pattern-analysis R script.
- Execute multistep pattern analysis pipeline, generating:
  1) Procrustes analysis between protein/phosphoprotein spaces and sample groups.
  2) Molecular subtyping and clustering (NMF-based subtype discovery).
  3) WGCNA co-expression module detection for both protein and phosphoprotein data.
- Perform comprehensive input validation (required paths, R script existence, metadata handling).
- Sanity-check output directories after execution to confirm results were generated.
- Provide debugging utilities (`get_command_string()`) for command inspection.

Simplified execution model:
- Direct `Rscript` execution (no conda environment activation).
- Assumes `Rscript` is available on the system PATH across all platforms.
- Cross-platform path handling via `to_r_path()` utility.

Typical usage:
- Construct `PatternAnalysisConfig` with required paths/parameters.
- Create `PatternAnalysis(config)` and call `run()` to execute the complete
  pattern analysis workflow (steps 1.1 → 1.2 → 1.3).
"""


from ..config import PatternAnalysisConfig
import pandas as pd
import subprocess
from pathlib import Path


def to_r_path(path: str) -> str:
    """Convert a local path to an R-friendly path string."""
    return str(Path(path)).replace('\\', '/')

# PatternAnalysis
class PatternAnalysis:
    """
    Pattern analysis wrapper.

    Responsibilities
    - Build an `Rscript` command from `PatternAnalysisConfig`.
    - Execute the R workflow to generate:
      1) Procrustes analysis
      2) Molecular subtyping (e.g., NMF)
      3) WGCNA co-expression modules
    """
    def __init__(self, config: PatternAnalysisConfig):
        self.config = config

    def build_command(self):
        """Build the `Rscript` command for the pattern analysis workflow."""
        cmd = [
            "Rscript",
            str(self.config.script_path / Path(self.config.pattern_analysis_rscript_path)),
            "--script_path", to_r_path(self.config.script_path),
            "--profile", to_r_path(self.config.profile),
            "--phosfile", to_r_path(self.config.phosfile),
            "--metadatafile", to_r_path(self.config.metadatafile),
            "--sample_group", to_r_path(self.config.sample_group),
            "--outdir", to_r_path(self.config.outdir),
            "--omics1_name", str(self.config.omics1_name),
            "--omics2_name", str(self.config.omics2_name),
            "--dim_rd_mtd", str(self.config.dim_rd_mtd),
            "--NMF_pro_filter_num", str(self.config.NMF_pro_filter_num),
            "--NMF_phos_filter_num", str(self.config.NMF_phos_filter_num),
            "--WGCNA_pro_filter_num", str(self.config.WGCNA_pro_filter_num),
            "--WGCNA_phos_filter_num", str(self.config.WGCNA_phos_filter_num),
            "--protein_SoftPower", str(self.config.protein_SoftPower),
            "--phosphoprotein_SoftPower", str(self.config.phosphoprotein_SoftPower),
            "--protein_RsquareCut", str(self.config.protein_RsquareCut_val),
            "--phosphoprotein_RsquareCut", str(self.config.phosphoprotein_RsquareCut_val),
            "--protein_cor_method", str(self.config.protein_cor_method),
            "--protein_corFun_tmp", str(self.config.protein_corFun_tmp),
            "--protein_cluster_method", str(self.config.protein_cluster_method),
            "--protein_corOptions", str(self.config.protein_corOptions_str),
            "--protein_networkType", str(self.config.protein_networkType),
            "--protein_mergingThresh", str(self.config.protein_mergingThresh),
            "--protein_minModuleSize", str(self.config.protein_minModuleSize),
            "--phosphoprotein_cor_method", str(self.config.phosphoprotein_cor_method),
            "--phosphoprotein_corFun_tmp", str(self.config.phosphoprotein_corFun_tmp),
            "--phosphoprotein_cluster_method", str(self.config.phosphoprotein_cluster_method),
            "--phosphoprotein_corOptions", str(self.config.phosphoprotein_corOptions_str),
            "--phosphoprotein_networkType", str(self.config.phosphoprotein_networkType),
            "--phosphoprotein_mergingThresh", str(self.config.phosphoprotein_mergingThresh),
            "--phosphoprotein_minModuleSize", str(self.config.phosphoprotein_minModuleSize),
            "--module_cor_threshhold", str(self.config.module_cor_threshhold),
            "--module_cor_p_adj", str(self.config.module_cor_p_adj),
            "--module_cor_method", str(self.config.module_cor_method),
            "--PA_group1_color", self.config.PA_group1_color,
            "--PA_group2_color", self.config.PA_group2_color,
            "--pro_ME_color", self.config.WGCNA_pro_ME_color,
            "--phos_ME_color", self.config.WGCNA_phos_ME_color,
            "--pro_color", self.config.WGCNA_pro_color,
            "--phos_color", self.config.WGCNA_phos_color,
            "--pheatmap_color", ";".join(self.config.WGCNA_pheatmap_color),
            # "--verbose"
        ]
        return cmd

    def validate_config(self):
        """Validate required configuration fields and input file existence."""
        required_attributes = [
            'script_path', 'profile', 'phosfile', 'sample_group',
            'outdir', 'omics1_name', 'omics2_name', 'pattern_analysis_rscript_path'
        ]
        
        for attr in required_attributes:
            if not hasattr(self.config, attr) or getattr(self.config, attr) is None:
                raise ValueError(f"Missing required configuration: {attr}")
        
        # Validate input files (paths) exist when they are provided as strings/Paths
        for file_attr in ['profile', 'phosfile', 'sample_group']:
            file_path = getattr(self.config, file_attr)
            if isinstance(file_path, (str, Path)) and not Path(file_path).exists():
                raise FileNotFoundError(f"File not found for {file_attr}: {file_path}")
        
        if isinstance(self.config.metadatafile, pd.DataFrame):
            has_metadata = not self.config.metadatafile.empty
        elif isinstance(self.config.metadatafile, (str, Path)):
            has_metadata = bool(str(self.config.metadatafile).strip())
        elif self.config.metadatafile is None:
            has_metadata = False
        else:
            has_metadata = False

        if not has_metadata:
            self.config.metadatafile = 'None'
        
        # Validate the R script exists
        script_full_path = self.config.script_path / Path(self.config.pattern_analysis_rscript_path)
        if not script_full_path.exists():
            raise FileNotFoundError(f"R script not found: {script_full_path}")

    def run(self):
        """Run the pattern analysis pipeline via the configured R script."""
        try:
            print("🔍 Starting pattern analysis...")
            
            # 1) Validate config
            self.validate_config()
            
            # 2) Build command
            cmd = self.build_command()
            
            # 3) Prepare execution command 
            full_cmd = cmd
            print(f"📋 Executing command:\n{' '.join(full_cmd)}")
            result = subprocess.run(full_cmd, capture_output=True, text=True)
            
            # 4) Handle results
            if result.returncode == 0:
                print("✅ Pattern analysis completed")
                if result.stdout:
                    print(f"📄 Output:\n{result.stdout}")
            else:
                print("❌ Pattern analysis failed")
                print(f"📄 STDOUT:\n{result.stdout}")
                print(f"❌ STDERR:\n{result.stderr}")
                return False
            
            # 5) Sanity-check output folders
            output_dirs = [
                self.config.outdir / "1.1Procrustes",
                self.config.outdir / "1.2Molecular_Subtyping", 
                self.config.outdir / "1.3WGCNA"
            ]
            
            for output_dir in output_dirs:
                if output_dir.exists() and any(output_dir.iterdir()):
                    print(f"✅ Output directory contains results: {output_dir}")
                else:
                    print(f"⚠️  Output directory is empty or missing: {output_dir}")
            
            return True
                
        except Exception as e:
            print(f"❌ Pattern analysis error: {e}")
            return False

    def get_command_string(self):
        """Return the assembled command string (useful for debugging)."""
        cmd = self.build_command()
        return ' '.join(cmd)