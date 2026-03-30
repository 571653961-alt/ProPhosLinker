
"""Analysis configuration definitions for ProPhosLinker.

Includes settings for pattern analysis, differential analysis, and functional analysis.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar, Union, List, Literal, Dict, Optional
import pandas as pd
from . import BasicConfig, ScriptConfig, Neo4jConfig, ResultDirConfig,ColorConfig

@dataclass
class PatternAnalysisConfig:
    """Pattern analysis configuration class."""
    
    # Input file paths (required when instantiating)
    profile: Union[str, Path, pd.DataFrame]
    phosfile: Union[str, Path, pd.DataFrame]
    sample_group: Union[str, Path, pd.DataFrame]
    metadatafile: Union[str, Path, pd.DataFrame] = None

    # Output directory
    outdir: Optional[Union[str, Path]] = None
    
    # Script path
    script_path: Optional[str] = None # ScriptConfig.script_path
    
    # Omics names
    omics1_name: Optional[str] = None #BasicConfig.omics1_name
    omics2_name: Optional[str] = None #BasicConfig.omics2_name
    
    # 1.1 PA parameters
    dim_rd_mtd: Literal["PCA", "PCoA"] = 'PCA'
    
    # 1.2 NMF parameters
    NMF_pro_filter_num: int = 3000
    NMF_phos_filter_num: int = 3000
    
    # 1.3 WGCNA parameters
    WGCNA_pro_filter_num: int = 5000
    WGCNA_phos_filter_num: int = 5000
    
    # Protein network parameters
    protein_cor_method: Literal["pearson", "spearman", "kendall"] = 'spearman'
    protein_corFun_tmp: Literal["bicor", "cor"] = 'bicor'
    protein_cluster_method: Literal["average", "complete", "ward.D"] = 'average'
    protein_corOptions_str: str = 'pairwise.complete.obs'   # cannot be changed, default
    protein_networkType: Literal["signed", "unsigned", "signed hybrid"] = 'signed'
    protein_RsquareCut_val: float = 0.89
    protein_mergingThresh: float = 0.2
    protein_minModuleSize: int = 10
    protein_SoftPower: int = 6
    
    # Phosphoprotein network parameters
    phosphoprotein_cor_method: Literal["pearson", "spearman", "kendall"] = 'spearman'
    phosphoprotein_corFun_tmp: Literal["bicor", "cor"] = 'bicor'
    phosphoprotein_cluster_method: Literal["average", "complete", "ward.D"] = 'average'
    phosphoprotein_corOptions_str: str = 'pairwise.complete.obs'
    phosphoprotein_networkType: Literal["signed", "unsigned", "signed hybrid"] = 'signed'
    phosphoprotein_RsquareCut_val: float = 0.89
    phosphoprotein_mergingThresh: float = 0.2
    phosphoprotein_minModuleSize: int = 10
    phosphoprotein_SoftPower: int = 6
    
    # Module-related parameters
    module_cor_threshhold: float = 0.5
    module_cor_p_adj: float = 0.05
    module_cor_method: Literal["pearson", "spearman", "kendall"] = 'spearman'
    
    # Output directory and script path
    pattern_analysis_rscript_path: Optional[Union[str, Path]] = None #ScriptConfig.get_full_path('pattern_analysis_rscript_path')
    
    # Plotting parameters
    # 1.1 PA colors
    PA_group1_color: str = ColorConfig.PA_group1_color
    PA_group2_color: str = ColorConfig.PA_group2_color
    
    # 1.3 WGCNA colors
    WGCNA_pro_ME_color: str = ColorConfig.WGCNA_pro_ME_color
    WGCNA_phos_ME_color: str = ColorConfig.WGCNA_phos_ME_color
    WGCNA_pro_color: str = ColorConfig.WGCNA_pro_color
    WGCNA_phos_color: str = ColorConfig.WGCNA_phos_color
    WGCNA_pheatmap_color: ClassVar[list] = ColorConfig.WGCNA_pheatmap_color

    def __post_init__(self):
        # outdir: if not provided, use BasicConfig + ResultDirConfig
        if self.outdir is None:
            base = getattr(BasicConfig, "outdir", None)
            if base is None:
                # fallback to cwd if BasicConfig not set
                self.outdir = Path.cwd()
            else:
                self.outdir = Path(base) / getattr(ResultDirConfig, "pattern_analysis", "pattern_analysis")
        else:
            self.outdir = Path(self.outdir)

        # script_path: prefer instance value, then ScriptConfig
        if not self.script_path:
            self.script_path = getattr(ScriptConfig, "script_path", "./scripts")

        # Omics names: prefer instance value, then BasicConfig
        if not self.omics1_name:
            self.omics1_name = getattr(BasicConfig, "omics1_name", "Pro")
        if not self.omics2_name:
            self.omics2_name = getattr(BasicConfig, "omics2_name", "Phos")

        if not self.pattern_analysis_rscript_path:
            self.pattern_analysis_rscript_path =  self.group_comparing = getattr(ScriptConfig, "pattern_analysis_rscript_path", "1.pattern_analysis.R")


    @classmethod
    def to_dict(cls) -> dict:
        """Convert configuration to a dictionary."""
        return {
            # Input files
            'profile': str(cls.profile) if isinstance(cls.profile, (str, Path)) else 'DataFrame',
            'phosfile': str(cls.phosfile) if isinstance(cls.phosfile, (str, Path)) else 'DataFrame',
            'sample_group': str(cls.sample_group) if isinstance(cls.sample_group, (str, Path)) else 'DataFrame',
            'metadatafile': str(cls.metadatafile) if isinstance(cls.metadatafile, (str, Path)) else 'DataFrame',
            
            # Output directory
            'outdir': str(cls.outdir),
            
            # Script path
            'script_path': cls.script_path,
            
            # Omics names
            'omics1_name': cls.omics1_name,
            'omics2_name': cls.omics2_name,
            
            # 1.1 PA parameters
            'dim_rd_mtd': cls.dim_rd_mtd,
            
            # 1.2 NMF parameters
            'NMF_pro_filter_num': cls.NMF_pro_filter_num,
            'NMF_phos_filter_num': cls.NMF_phos_filter_num,
            
            # 1.3 WGCNA parameters - protein
            'WGCNA_pro_filter_num': cls.WGCNA_pro_filter_num,
            'protein_cor_method': cls.protein_cor_method,
            'protein_corFun_tmp': cls.protein_corFun_tmp,
            'protein_cluster_method': cls.protein_cluster_method,
            'protein_corOptions_str': cls.protein_corOptions_str,
            'protein_networkType': cls.protein_networkType,
            'protein_RsquareCut_val': cls.protein_RsquareCut_val,
            'protein_mergingThresh': cls.protein_mergingThresh,
            'protein_minModuleSize': cls.protein_minModuleSize,
            'protein_SoftPower': cls.protein_SoftPower,
            
            # 1.3 WGCNA parameters - phosphoprotein
            'WGCNA_phos_filter_num': cls.WGCNA_phos_filter_num,
            'phosphoprotein_cor_method': cls.phosphoprotein_cor_method,
            'phosphoprotein_corFun_tmp': cls.phosphoprotein_corFun_tmp,
            'phosphoprotein_cluster_method': cls.phosphoprotein_cluster_method,
            'phosphoprotein_corOptions_str': cls.phosphoprotein_corOptions_str,
            'phosphoprotein_networkType': cls.phosphoprotein_networkType,
            'phosphoprotein_RsquareCut_val': cls.phosphoprotein_RsquareCut_val,
            'phosphoprotein_mergingThresh': cls.phosphoprotein_mergingThresh,
            'phosphoprotein_minModuleSize': cls.phosphoprotein_minModuleSize,
            'phosphoprotein_SoftPower': cls.phosphoprotein_SoftPower,
            
            # Module-related parameters
            'module_cor_threshhold': cls.module_cor_threshhold,
            'module_cor_p_adj': cls.module_cor_p_adj,
            'module_cor_method': cls.module_cor_method,
            
            # Output directory and script path
            # 'oudir_default_name': str(cls.oudir_default_name),
            'pattern_analysis_rscript_path': str(cls.pattern_analysis_rscript_path),
            
            # Plotting parameters
            # 1.1 PA colors
            'PA_group1_color': cls.PA_group1_color,
            'PA_phos_color': cls.PA_phos_color,
            
            # 1.3 WGCNA colors
            'WGCNA_pro_ME_color': cls.WGCNA_pro_ME_color,
            'WGCNA_phos_ME_color': cls.WGCNA_phos_ME_color,
            'WGCNA_pro_color': cls.WGCNA_pro_color,
            'WGCNA_phos_color': cls.WGCNA_phos_color,
            'WGCNA_pheatmap_color': cls.WGCNA_pheatmap_color
        }
    
    @classmethod
    def update_config(cls, **kwargs):
        """Update configuration parameters."""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
    
    @classmethod
    def get_analysis_sections(cls) -> dict:
        """Get analysis module groupings."""
        return {
            'pattern_analysis_1_1': [
                'dim_rd_mtd','PA_group1_color', 'PA_phos_color'
            ],
            'nmf_analysis_1_2': [
                'NMF_pro_filter_num', 'NMF_phos_filter_num'
            ],
            'wgcna_analysis_1_3': [
                'WGCNA_pro_filter_num', 'WGCNA_phos_filter_num',
                'protein_cor_method', 'protein_corFun_tmp', 'protein_cluster_method',
                'protein_corOptions_str', 'protein_networkType', 'protein_RsquareCut_val',
                'protein_mergingThresh', 'protein_minModuleSize', 'protein_SoftPower',
                'phosphoprotein_cor_method', 'phosphoprotein_corFun_tmp', 'phosphoprotein_cluster_method',
                'phosphoprotein_corOptions_str', 'phosphoprotein_networkType', 'phosphoprotein_RsquareCut_val',
                'phosphoprotein_mergingThresh', 'phosphoprotein_minModuleSize', 'phosphoprotein_SoftPower',
                'module_cor_threshhold', 'module_cor_p_adj', 'module_cor_method',
                'WGCNA_pro_ME_color', 'WGCNA_phos_ME_color',
                'WGCNA_pro_color', 'WGCNA_phos_color', 'WGCNA_pheatmap_color'
            ]
        }
    
    @classmethod
    def validate_input_files(cls) -> bool:
        """Validate that required input files exist."""
        import os
        
        files_to_check = [
            cls.profile,
            cls.phosfile, 
            cls.sample_group,
            cls.metadatafile
        ]
        
        missing_files = []
        for file_path in files_to_check:
            if isinstance(file_path, (str, Path)) and not os.path.exists(file_path):
                missing_files.append(str(file_path))
        
        if missing_files:
            print("❌ Missing input files:")
            for file in missing_files:
                print(f"   - {file}")
            return False
        
        print("✅ All input files exist")
        return True
    
    @classmethod
    def get_wgcna_parameters(cls, omics_type: Literal["protein", "phosphoprotein"]) -> dict:
        """Get WGCNA parameters for a specified omics type."""
        prefix = omics_type.lower()
        return {
            'filter_num': getattr(cls, f'WGCNA_{prefix}_filter_num'),
            'cor_method': getattr(cls, f'{prefix}_cor_method'),
            'corFun_tmp': getattr(cls, f'{prefix}_corFun_tmp'),
            'cluster_method': getattr(cls, f'{prefix}_cluster_method'),
            'corOptions_str': getattr(cls, f'{prefix}_corOptions_str'),
            'networkType': getattr(cls, f'{prefix}_networkType'),
            'RsquareCut_val': getattr(cls, f'{prefix}_RsquareCut_val'),
            'mergingThresh': getattr(cls, f'{prefix}_mergingThresh'),
            'minModuleSize': getattr(cls, f'{prefix}_minModuleSize'),
            'SoftPower': getattr(cls, f'{prefix}_SoftPower')
        }

"""
Differential analysis configuration settings for the project.
"""

@dataclass
class DifferentialConfig:
    """Differential analysis configuration class."""
    
    # Input file paths (required when instantiating)
    profile: Union[str, Path, pd.DataFrame]
    phosfile: Union[str, Path, pd.DataFrame]
    mappingfile: Union[str, Path, pd.DataFrame]
    sample_group: Union[str, Path, pd.DataFrame]
    
    # Optional input files
    # NOTE: Do not use ClassVar — dataclass will exclude ClassVar fields from __init__,
    # which would cause DifferentialConfig(prodiff=...) to raise "unexpected keyword argument".
    prodiff: Optional[Union[str, Path, pd.DataFrame]] = None
    phosdiff: Optional[Union[str, Path, pd.DataFrame]] = None
    
    # Output directory
    outdir: Optional[Union[str, Path]] = None
    # oudir_default_name: str = ResultDirConfig.differential_analysis
    
    # Script path
    script_path: str = ScriptConfig.script_path
    
    # Omics names
    omics1_name: str = None
    omics2_name: str = None
    
    # Identification type
    identified_type: Literal['UNIPROT', 'SYMBOL'] = None
    
    # Group comparison
    group_comparing: str = None
    
    # 2.1 Differential analysis parameters
    pro_log2FC: float = BasicConfig.pro_log2FC
    phos_log2FC: float = BasicConfig.phos_log2FC
    diff_p_val: float = BasicConfig.pro_diff_q_val
    omics_differential_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.get_full_path('omics_differential_rscript_path') #ScriptConfig.omics_differential_rscript_path
    
    # 2.1 Plotting parameters
    omics_diff_quadrant_pro_color: str = ColorConfig.omics_diff_quadrant_pro_color
    omics_diff_quadrant_phos_color: str = ColorConfig.omics_diff_quadrant_phos_color
    
    # 2.2 Phosphorylation rate subtyping parameters
    phos_rate_FC: float = BasicConfig.phosRate_log2FC
    quantile_rank: int = 100
    quantile_rank_step: int = 1
    frac_param: float = 1
    cluster_num: int = 9
    plot_sites_num: int = 7   ## This seems unused
    plot_site_ids: ClassVar[List] = None   ## This seems unused
    vs_total_phos_rate_upper: float = 0.6
    vs_total_phos_rate_lower: float = 0.4
    top_tail_diff_site_num: int = 5
    
    # 2.2 Plotting parameters
    raw_pro_point_color: str = ColorConfig.raw_pro_point_color
    fitted_pro_point_color: str = ColorConfig.fitted_pro_point_color
    fitted_pro_line_color: str = ColorConfig.fitted_pro_line_color
    pro_fc_line_color: str = ColorConfig.pro_fc_line_color
    total_phos_rate_line_color: str = ColorConfig.total_phos_rate_line_color
    color_lst: ClassVar[List[str]] = ColorConfig.color_lst
    Mfuzz_pluscolor: str = ColorConfig.Mfuzz_pluscolor
    Mfuzz_minuscolor: str = ColorConfig.Mfuzz_minuscolor
    Mfuzz_allcolor: str = ColorConfig.Mfuzz_allcolor
    Mfuzz_onlyupcolor: str = ColorConfig.Mfuzz_onlyupcolor
    Mfuzz_onlydowncolor: str = ColorConfig.Mfuzz_onlydowncolor
    Mfuzz_updowncolor: str = ColorConfig.Mfuzz_updowncolor
    Mfuzz_notsigcolor: str = ColorConfig.Mfuzz_notsigcolor
    phosRate_quantile_subtyping_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.get_full_path('phosRate_quantile_subtyping_rscript_path') #ScriptConfig.phosRate_quantile_subtyping_rscript_path
    
    # 2.3 Differential network parameters
    differential_network_diff_pro_path: Union[str, Path] = None
    differential_network_diff_phos_path: Union[str, Path] = None
    differential_network_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.get_full_path('differential_network_rscript_path') #ScriptConfig.differential_network_rscript_path
    differential_network_filter_num: int = 3000
    differential_network_FC_threshold: float = 2
    differential_network_p_threshold: float = 0.05
    differential_network_p_value_type: Literal["q_value", "P_value"] = "q_value"
    differential_network_nBoots: int = 20
    differential_network_bootnet_R_threshold: float = 0.3
    differential_network_nCores: int = 6
    differential_network_stability_threshold: float = 0.4
    differential_network_cor_method: Literal["pearson", "spearman", "kendall"] = "spearman"
    differential_network_edge_FC_threshold: float = 2
    differential_network_edge_p_threshold: float = 0.05
    differential_network_max_subnet_num: int = 6
    differential_network_R_threshold: float = 0.5
    
    # 2.3 Differential network plotting parameters
    differential_network_edge_color_pos: str = ColorConfig.differential_network_edge_color_pos
    differential_network_edge_color_neg: str = ColorConfig.differential_network_edge_color_neg
    differential_network_Enhanced_in_N: str = ColorConfig.differential_network_Enhanced_in_N
    differential_network_Enhanced_in_T: str = ColorConfig.differential_network_Enhanced_in_T
    differential_network_Only_in_N: str = ColorConfig.differential_network_Only_in_N
    differential_network_Only_in_T: str = ColorConfig.differential_network_Only_in_T
    differential_network_Conflict_relation: str = ColorConfig.differential_network_Conflict_relation
    differential_network_gradientn_color: ClassVar[List[str]] = ColorConfig.differential_network_gradientn_color
    differential_network_function_enrichment_color_gradient_low: str = ColorConfig.omics_function_enrichment_color_gradient_low
    differential_network_function_enrichment_color_gradient_high: str = ColorConfig.omics_function_enrichment_color_gradient_high

    def __post_init__(self):
        # outdir: if not provided, use BasicConfig + ResultDirConfig
        if self.outdir is None:
            base = getattr(BasicConfig, "outdir", None)
            if base is None:
                # fallback to cwd if BasicConfig not set
                self.outdir = Path.cwd()
            else:
                self.outdir = Path(base) / getattr(ResultDirConfig, "differential_analysis", "differential_analysis")
        else:
            self.outdir = Path(self.outdir)

        # script_path: prefer instance value, then ScriptConfig
        if not self.script_path:
            self.script_path = getattr(ScriptConfig, "script_path", "./scripts")

        # omics names: prefer instance value, then BasicConfig
        if not self.omics1_name:
            self.omics1_name = getattr(BasicConfig, "omics1_name", "Pro")
        if not self.omics2_name:
            self.omics2_name = getattr(BasicConfig, "omics2_name", "Phos")
        if not self.group_comparing:
            self.group_comparing = BasicConfig.group_comparing
        if not self.identified_type:
            self.identified_type = BasicConfig.identified_type

        if not self.omics_differential_rscript_path:
            self.omics_differential_rscript_path =  getattr(ScriptConfig, "omics_differential_rscript_path", "2.1omics_differential_analysis.R")

    @classmethod
    def to_dict(cls) -> dict:
        """Convert configuration to a dictionary."""
        return {
            # Input files
            'profile': str(cls.profile) if isinstance(cls.profile, (str, Path)) else 'DataFrame',
            'phosfile': str(cls.phosfile) if isinstance(cls.phosfile, (str, Path)) else 'DataFrame',
            'mappingfile': str(cls.mappingfile) if isinstance(cls.mappingfile, (str, Path)) else 'DataFrame',
            'sample_group': str(cls.sample_group) if isinstance(cls.sample_group, (str, Path)) else 'DataFrame',
            'prodiff': str(cls.prodiff) if isinstance(cls.prodiff, (str, Path)) else 'DataFrame',
            'phosdiff': str(cls.phosdiff) if isinstance(cls.phosdiff, (str, Path)) else 'DataFrame',
            
            # Output directory
            'outdir': str(cls.outdir),
            'oudir_default_name': cls.oudir_default_name,
            
            # Script path
            'script_path': cls.script_path,
            
            # Omics names
            'omics1_name': cls.omics1_name,
            'omics2_name': cls.omics2_name,
            
            # Identification type
            'identified_type': cls.identified_type,
            
            # Group comparison
            'group_comparing': cls.group_comparing,
            
            # 2.1 Differential analysis parameters
            'pro_log2FC': cls.pro_log2FC,
            'phos_log2FC': cls.phos_log2FC,
            'diff_p_val': cls.diff_p_val,
            'omics_differential_rscript_path': str(cls.omics_differential_rscript_path),
            
            # 2.1 Plotting parameters
            'omics_diff_quadrant_pro_color': cls.omics_diff_quadrant_pro_color,
            'omics_diff_quadrant_phos_color': cls.omics_diff_quadrant_phos_color,
            
            # 2.2 Phosphorylation rate subtyping parameters
            'phos_rate_FC': cls.phos_rate_FC,
            'quantile_rank': cls.quantile_rank,
            'quantile_rank_step': cls.quantile_rank_step,
            'frac_param': cls.frac_param,
            'cluster_num': cls.cluster_num,
            'plot_sites_num': cls.plot_sites_num,
            'plot_site_ids': cls.plot_site_ids,
            'vs_total_phos_rate_upper': cls.vs_total_phos_rate_upper,
            'vs_total_phos_rate_lower': cls.vs_total_phos_rate_lower,
            'top_tail_diff_site_num': cls.top_tail_diff_site_num,
            
            # 2.2 Plotting parameters
            'raw_pro_point_color': cls.raw_pro_point_color,
            'fitted_pro_point_color': cls.fitted_pro_point_color,
            'fitted_pro_line_color': cls.fitted_pro_line_color,
            'pro_fc_line_color': cls.pro_fc_line_color,
            'total_phos_rate_line_color': cls.total_phos_rate_line_color,
            'color_lst': cls.color_lst,
            'Mfuzz_pluscolor': cls.Mfuzz_pluscolor,
            'Mfuzz_minuscolor': cls.Mfuzz_minuscolor,
            'Mfuzz_allcolor': cls.Mfuzz_allcolor,
            'Mfuzz_onlyupcolor': cls.Mfuzz_onlyupcolor,
            'Mfuzz_onlydowncolor': cls.Mfuzz_onlydowncolor,
            'Mfuzz_updowncolor': cls.Mfuzz_updowncolor,
            'Mfuzz_notsigcolor': cls.Mfuzz_notsigcolor,
            'phosRate_quantile_subtyping_rscript_path': str(cls.phosRate_quantile_subtyping_rscript_path),
            
            # 2.3 Differential network parameters
            'differential_network_rscript_path': str(cls.differential_network_rscript_path),
            'differential_network_filter_num': cls.differential_network_filter_num,
            'differential_network_FC_threshold': cls.differential_network_FC_threshold,
            'differential_network_p_threshold': cls.differential_network_p_threshold,
            'differential_network_p_value_type': cls.differential_network_p_value_type,
            'differential_network_nBoots': cls.differential_network_nBoots,
            'differential_network_bootnet_R_threshold': cls.differential_network_bootnet_R_threshold,
            'differential_network_nCores': cls.differential_network_nCores,
            'differential_network_stability_threshold': cls.differential_network_stability_threshold,
            'differential_network_cor_method': cls.differential_network_cor_method,
            'differential_network_edge_FC_threshold': cls.differential_network_edge_FC_threshold,
            'differential_network_edge_p_threshold': cls.differential_network_edge_p_threshold,
            'differential_network_max_subnet_num': cls.differential_network_max_subnet_num,
            'differential_network_R_threshold': cls.differential_network_R_threshold,
            
            # 2.3 Differential network plotting parameters
            'differential_network_edge_color_pos': cls.differential_network_edge_color_pos,
            'differential_network_edge_color_neg': cls.differential_network_edge_color_neg,
            'differential_network_Enhanced_in_N': cls.differential_network_Enhanced_in_N,
            'differential_network_Enhanced_in_T': cls.differential_network_Enhanced_in_T,
            'differential_network_Only_in_N': cls.differential_network_Only_in_N,
            'differential_network_Only_in_T': cls.differential_network_Only_in_T,
            'differential_network_Conflict_relation': cls.differential_network_Conflict_relation,
            'differential_network_gradientn_color': cls.differential_network_gradientn_color,
            'differential_network_function_enrichment_color_gradient_low': cls.differential_network_function_enrichment_color_gradient_low,
            'differential_network_function_enrichment_color_gradient_high': cls.differential_network_function_enrichment_color_gradient_high,
        }
    
    @classmethod
    def update_config(cls, **kwargs):
        """Update configuration parameters."""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
    
    @classmethod
    def get_analysis_sections(cls) -> dict:
        """Get analysis module groupings."""
        return {
            'differential_analysis_2_1': [
                'pro_log2FC', 'phos_log2FC', 'diff_p_val',
                'omics_diff_quadrant_pro_color', 'omics_diff_quadrant_phos_color'
            ],
            'phosrate_subtyping_2_2': [
                'phos_rate_FC', 'quantile_rank', 'quantile_rank_step', 'frac_param',
                'cluster_num', 'plot_sites_num', 'plot_site_ids',
                'vs_total_phos_rate_upper', 'vs_total_phos_rate_lower', 'top_tail_diff_site_num',
                'raw_pro_point_color', 'fitted_pro_point_color', 'fitted_pro_line_color',
                'pro_fc_line_color', 'total_phos_rate_line_color', 'color_lst',
                'Mfuzz_pluscolor', 'Mfuzz_minuscolor', 'Mfuzz_allcolor',
                'Mfuzz_onlyupcolor', 'Mfuzz_onlydowncolor', 'Mfuzz_updowncolor', 'Mfuzz_notsigcolor'
            ],
            'differential_network_2_3': [
                'differential_network_filter_num', 'differential_network_FC_threshold',
                'differential_network_p_threshold', 'differential_network_p_value_type',
                'differential_network_nBoots', 'differential_network_bootnet_R_threshold',
                'differential_network_nCores', 'differential_network_stability_threshold',
                'differential_network_cor_method', 'differential_network_edge_FC_threshold',
                'differential_network_edge_p_threshold', 'differential_network_max_subnet_num',
                'differential_network_R_threshold', 'differential_network_edge_color_pos',
                'differential_network_edge_color_neg', 'differential_network_Enhanced_in_N',
                'differential_network_Enhanced_in_T', 'differential_network_Only_in_N',
                'differential_network_Only_in_T', 'differential_network_Conflict_relation',
                'differential_network_gradientn_color',
                'differential_network_function_enrichment_color_gradient_low','differential_network_function_enrichment_color_gradient_high'
            ]
        }


"""
Functional analysis configuration settings for the project.
"""
@dataclass
class FunctionalConfig:
    """Functional analysis configuration class."""
    
    # Input file paths (required when instantiating)
    profile: Union[str, Path, pd.DataFrame]
    phosfile: Union[str, Path, pd.DataFrame]
    mappingfile: Union[str, Path, pd.DataFrame]
    sample_group: Union[str, Path, pd.DataFrame]
    
    # Optional input files
    prodiff: Union[str, Path, pd.DataFrame]
    phosdiff: Union[str, Path, pd.DataFrame]
    
    # Output directory
    outdir: ClassVar[Union[str, Path]] = None
    
    # Script path
    script_path: str = None
    
    # Omics names
    omics1_name: str = None
    omics2_name: str = None
    
    # Identification type
    identified_type: Literal['UNIPROT', 'SYMBOL'] = None
    
    # Group comparison
    group_comparing: str = None
    
    # 3.1 Functional enrichment analysis parameters
    network_log2FC: float = None
    GO_showCategory: int = 6
    KEGG_showCategory: int = 18
    diff_p_adj: float = 0.05
    pvalueCutoff: float = 0.05
    
    omics_function_enrichment_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.omics_function_enrichment_rscript_path
    
    # 3.1 Functional enrichment plotting parameters
    omics_function_enrichment_color_gradient_low: str = ColorConfig.omics_function_enrichment_color_gradient_low
    omics_function_enrichment_color_gradient_high: str = ColorConfig.omics_function_enrichment_color_gradient_high
    
    # 3.2 Knowledge graph functional network parameters
    uri: str = Neo4jConfig.uri
    username: str = Neo4jConfig.username
    password: str = Neo4jConfig.password
    kgbased_functional_network_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.kgbased_functional_network_rscript_path
    kgbased_functional_network_community_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.kgbased_functional_network_community_rscript_path
    kgbased_network_disease: str = BasicConfig.disease
    kgbased_network_top_n: int = 6
    kgbased_network_analysis_mode : Literal["ALL",    # no public interface; kept in utilities
        "KS","LR","TF",
        "DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'DEP'
    kgbased_network_top_nodes_visualization_num: int = 20
    kgbased_network_max_phosphoSite_displayed: int = 5
    kgbased_network_visualization_mode : Literal["ALL","KS","LR","TF","DEP", "DEPUP", "DEPDOWN", "SEEDNODE"] = 'ALL'
    kgbased_network_node_filtering : Literal["betweenness", "degree",  "closeness", "harmonic","eigenvector", "pagerank", "alpha", "hub",  "authority"] = 'betweenness'
    kgbased_network_layout : Literal["fr","kk","dh","stress","tree","gem","graphopt", "lgl","circle","grid"] = 'kk'
    kgbased_network_comm_detection : Literal["louvain","leiden","infomap","fastgreedy","walktrap", "LPA"] ='infomap'
    
    # 3.2 Knowledge graph functional network plotting parameters
    kgbased_functional_network_node_up: str = ColorConfig.kgbased_functional_network_node_up
    kgbased_functional_network_node_down: str = ColorConfig.kgbased_functional_network_node_down
    kgbased_functional_network_node_nonsig: str = ColorConfig.kgbased_functional_network_node_nonsig
    kgbased_functional_network_node_notdet: str = ColorConfig.kgbased_functional_network_node_notdet
    kgbased_functional_network_node_disease_border_color: str = ColorConfig.kgbased_functional_network_node_disease_border_color
    kgbased_functional_network_node_notdisease_border_color: str = ColorConfig.kgbased_functional_network_node_notdisease_border_color
    kgbased_functional_network_function_enrich_color_gradient_low: str = ColorConfig.kgbased_functional_network_function_enrich_color_gradient_low
    kgbased_functional_network_function_enrich_color_gradient_high: str = ColorConfig.kgbased_functional_network_function_enrich_color_gradient_high
    kgbased_functional_network_edge_type_colors: ClassVar[Dict[str, str]] = ColorConfig.kgbased_functional_network_edge_type_colors
    
    # 3.3 Hub protein functional network parameters
    hub_protein_top_n: int = 6
    hub_protein_d_num: int = 5
    hub_protein_network_mode: str = "ALL"
    hub_protein_network_layout: Literal["fr", "kk", "dh", "stress", "tree", "gem", "graphopt", "lgl", "circle", "grid"] = 'kk'
    function_enrich_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.function_enrich_rscript_path
    hub_protein_functional_network_rscript_path: ClassVar[Union[str, Path]] = ScriptConfig.hub_protein_functional_network_rscript_path
    
    # 3.3 Hub protein functional network plotting parameters (uses the same color mapping as the knowledge graph)
    hub_protein_network_node_up: str = ColorConfig.hub_protein_network_node_up
    hub_protein_network_node_down: str = ColorConfig.hub_protein_network_node_down
    hub_protein_network_node_nonsig: str = ColorConfig.hub_protein_network_node_nonsig
    hub_protein_network_node_notdet: str = ColorConfig.hub_protein_network_node_notdet
    hub_protein_network_node_disease_border_color: str = ColorConfig.hub_protein_network_node_disease_border_color
    hub_protein_network_node_notdisease_border_color: str = ColorConfig.hub_protein_network_node_notdisease_border_color
    hub_protein_network_function_enrich_color_gradient_low: str = ColorConfig.hub_protein_network_function_enrich_color_gradient_low
    hub_protein_network_function_enrich_color_gradient_high: str = ColorConfig.hub_protein_network_function_enrich_color_gradient_high
    hub_protein_network_edge_type_colors: ClassVar[Dict[str, str]] = ColorConfig.hub_protein_network_edge_type_colors

    def __post_init__(self):
        # outdir: if not provided, use BasicConfig + ResultDirConfig
        if self.outdir is None:
            base = getattr(BasicConfig, "outdir", None)
            if base is None:
                # fallback to cwd if BasicConfig not set
                self.outdir = Path.cwd()
            else:
                self.outdir = Path(base) / getattr(ResultDirConfig, "functional_analysis", "functional_analysis")
        else:
            self.outdir = Path(self.outdir)

        # script_path: prefer instance value, then ScriptConfig
        if not self.script_path:
            self.script_path = getattr(ScriptConfig, "script_path", "./scripts")

        # omics names: prefer instance value, then BasicConfig
        if not self.omics1_name:
            self.omics1_name = getattr(BasicConfig, "omics1_name", "Pro")
        if not self.omics2_name:
            self.omics2_name = getattr(BasicConfig, "omics2_name", "Phos")
        if not self.group_comparing:
            self.group_comparing = BasicConfig.group_comparing
        if not self.identified_type:
            self.identified_type = BasicConfig.identified_type

        if not self.network_log2FC:
            self.network_log2FC = BasicConfig.network_log2FC
        if not self.username:
            self.username = Neo4jConfig.username
        if not self.password:
            self.self.password = Neo4jConfig.password
        if not self.kgbased_network_disease:
            self.kgbased_network_disease: str = BasicConfig.disease


    @classmethod
    def to_dict(cls) -> dict:
        """Convert configuration to a dictionary."""
        return {
            # Input files
            'profile': str(cls.profile) if isinstance(cls.profile, (str, Path)) else 'DataFrame',
            'phosfile': str(cls.phosfile) if isinstance(cls.phosfile, (str, Path)) else 'DataFrame',
            'mappingfile': str(cls.mappingfile) if isinstance(cls.mappingfile, (str, Path)) else 'DataFrame',
            'sample_group': str(cls.sample_group) if isinstance(cls.sample_group, (str, Path)) else 'DataFrame',
            'prodiff': str(cls.prodiff) if isinstance(cls.prodiff, (str, Path)) else 'DataFrame',
            'phosdiff': str(cls.phosdiff) if isinstance(cls.phosdiff, (str, Path)) else 'DataFrame',
            
            # Output directory
            'outdir': str(cls.outdir),
            'oudir_default_name': cls.oudir_default_name,
            
            # Script path
            'script_path': cls.script_path,
            
            # Omics names
            'omics1_name': cls.omics1_name,
            'omics2_name': cls.omics2_name,
            
            # Identification type
            'identified_type': cls.identified_type,
            
            # Group comparison
            'group_comparing': cls.group_comparing,
            
            # 3.1 Functional enrichment analysis parameters
            'network_log2FC': cls.network_log2FC,
            'diff_p_adj': cls.diff_p_adj,
            'GO_showCategory': cls.GO_showCategory,
            'KEGG_showCategory': cls.KEGG_showCategory,
            'pvalueCutoff': cls.pvalueCutoff,
            'omics_function_enrichment_rscript_path': str(cls.omics_function_enrichment_rscript_path),
            
            # 3.1 Functional enrichment plotting parameters
            'omics_function_enrichment_color_gradient_low': cls.omics_function_enrichment_color_gradient_low,
            'omics_function_enrichment_color_gradient_high': cls.omics_function_enrichment_color_gradient_high,
            
            # 3.2 Knowledge graph functional network parameters
            'uri': cls.uri,
            'username': cls.username,
            'password': '***' if cls.password else None,  # hide password
            'kgbased_functional_network_rscript_path': str(cls.kgbased_functional_network_rscript_path),
            'kgbased_functional_network_community_rscript_path': str(cls.kgbased_functional_network_community_rscript_path),
            'kgbased_network_disease': cls.kgbased_network_disease,
            'kgbased_network_top_n': cls.kgbased_network_top_n,
            
            # 3.2 Knowledge graph functional network plotting parameters
            'kgbased_functional_network_node_up': cls.kgbased_functional_network_node_up,
            'kgbased_functional_network_node_down': cls.kgbased_functional_network_node_down,
            'kgbased_functional_network_node_nonsig': cls.kgbased_functional_network_node_nonsig,
            'kgbased_functional_network_node_notdet': cls.kgbased_functional_network_node_notdet,
            'kgbased_functional_network_node_disease_border_color': cls.kgbased_functional_network_node_disease_border_color,
            'kgbased_functional_network_node_notdisease_border_color': cls.kgbased_functional_network_node_notdisease_border_color,
            'kgbased_functional_network_function_enrich_color_gradient_low': cls.kgbased_functional_network_function_enrich_color_gradient_low,
            'kgbased_functional_network_function_enrich_color_gradient_high': cls.kgbased_functional_network_function_enrich_color_gradient_high,
            'kgbased_functional_network_edge_type_colors': cls.kgbased_functional_network_edge_type_colors,
            
            # 3.3 Hub protein functional network parameters
            'hub_protein_top_n': cls.hub_protein_top_n,
            'hub_protein_d_num': cls.hub_protein_d_num,
            'hub_protein_network_mode': cls.hub_protein_network_mode,
            'hub_protein_network_layout': cls.hub_protein_network_layout,
            'function_enrich_rscript_path': str(cls.function_enrich_rscript_path),
            'hub_protein_functional_network_rscript_path': str(cls.hub_protein_functional_network_rscript_path),
            
            # 3.3 Hub protein functional network plotting parameters
            'hub_protein_network_node_up': cls.hub_protein_network_node_up,
            'hub_protein_network_node_down': cls.hub_protein_network_node_down,
            'hub_protein_network_node_nonsig': cls.hub_protein_network_node_nonsig,
            'hub_protein_network_node_notdet': cls.hub_protein_network_node_notdet,
            'hub_protein_network_node_disease_border_color': cls.hub_protein_network_node_disease_border_color,
            'hub_protein_network_node_notdisease_border_color': cls.hub_protein_network_node_notdisease_border_color,
            'hub_protein_network_function_enrich_color_gradient_low': cls.hub_protein_network_function_enrich_color_gradient_low,
            'hub_protein_network_function_enrich_color_gradient_high': cls.hub_protein_network_function_enrich_color_gradient_high,
            'hub_protein_network_edge_type_colors': cls.hub_protein_network_edge_type_colors
        }
    
    @classmethod
    def update_config(cls, **kwargs):
        """Update configuration parameters."""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
    
    @classmethod
    def get_analysis_sections(cls) -> dict:
        """Get analysis module groupings."""
        return {
            'functional_enrichment_3_1': [
                'network_log2FC', 'diff_p_adj', 'GO_showCategory', 'KEGG_showCategory', 'pvalueCutoff',
                'omics_function_enrichment_color_gradient_low', 'omics_function_enrichment_color_gradient_high'
            ],
            'kgbased_network_3_2': [
                'uri', 'username', 'password', 'kgbased_network_disease', 'kgbased_network_top_n',
                'kgbased_functional_network_node_up', 'kgbased_functional_network_node_down',
                'kgbased_functional_network_node_nonsig', 'kgbased_functional_network_node_notdet',
                'kgbased_functional_network_node_disease_border_color', 'kgbased_functional_network_node_notdisease_border_color',
                'kgbased_functional_network_function_enrich_color_gradient_low', 'kgbased_functional_network_function_enrich_color_gradient_high',
                'kgbased_functional_network_edge_type_colors'
            ],
            'hub_protein_network_3_3': [
                'hub_protein_top_n', 'hub_protein_d_num', 'hub_protein_network_mode', 'hub_protein_network_layout',
                'hub_protein_network_node_up', 'hub_protein_network_node_down',
                'hub_protein_network_node_nonsig', 'hub_protein_network_node_notdet',
                'hub_protein_network_node_disease_border_color', 'hub_protein_network_node_notdisease_border_color',
                'hub_protein_network_function_enrich_color_gradient_low', 'hub_protein_network_function_enrich_color_gradient_high',
                'hub_protein_network_edge_type_colors'
            ]
        }
    
    @classmethod
    def validate_input_files(cls) -> bool:
        """Validate that required input files exist."""
        import os
        
        files_to_check = [
            cls.profile,
            cls.phosfile,
            cls.mappingfile,
            cls.sample_group
        ]
        
        missing_files = []
        for file_path in files_to_check:
            if isinstance(file_path, (str, Path)) and not os.path.exists(file_path):
                missing_files.append(str(file_path))
        
        if missing_files:
            print("❌ Missing input files:")
            for file in missing_files:
                print(f"   - {file}")
            return False
        
        print("✅ All input files exist")
        return True
    
    @classmethod
    def get_neo4j_config(cls) -> dict:
        """Get Neo4j configuration."""
        return {
            'uri': cls.uri,
            'username': cls.username,
            'password': cls.password
        }
    
    @classmethod
    def get_network_colors(cls, network_type: Literal["kgbased", "hub_protein"]) -> dict:
        """Get network color configuration."""
        prefix = "kgbased_functional_network" if network_type == "kgbased" else "hub_protein_network"
        
        return {
            'node_up': getattr(cls, f'{prefix}_node_up'),
            'node_down': getattr(cls, f'{prefix}_node_down'),
            'node_nonsig': getattr(cls, f'{prefix}_node_nonsig'),
            'node_notdet': getattr(cls, f'{prefix}_node_notdet'),
            'node_disease_border_color': getattr(cls, f'{prefix}_node_disease_border_color'),
            'node_notdisease_border_color': getattr(cls, f'{prefix}_node_notdisease_border_color'),
            'function_enrich_color_gradient_low': getattr(cls, f'{prefix}_function_enrich_color_gradient_low'),
            'function_enrich_color_gradient_high': getattr(cls, f'{prefix}_function_enrich_color_gradient_high'),
            'edge_type_colors': getattr(cls, f'{prefix}_edge_type_colors')
        }