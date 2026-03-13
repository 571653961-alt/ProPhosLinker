#!/usr/bin/env python3
"""
DataPreprocessingConfig: Comprehensive Data Preprocessing Configuration

Dataclass configuration for the complete proteomics/phosphoproteomics preprocessing 
pipeline. Supports flexible input formats (file paths or DataFrames), intelligent 
default resolution from other config classes, and comprehensive validation of 
filtering/normalization/imputation parameters.

Key Features:
- Multi-format input support (str/Path/DataFrame)
- Auto-resolution of defaults from BasicConfig/ScriptConfig/ResultDirConfig
- 7 imputation methods (zero/min/mean/median/knn/gmm/MLE)
- 4 normalization methods (median/quantile/PQN/none)
- Missing value filtering thresholds
- Post-init validation and Path normalization
- Comprehensive to_dict() for downstream serialization

Pipeline Parameters:
1. Input files (5 required + 2 optional differential results)
2. Filtering (missing value ratios)
3. Imputation (method + KNN neighbors)
4. Normalization (omics-specific methods)
5. Output directory auto-management
"""


from dataclasses import dataclass
from typing import ClassVar, Union, Literal, Optional
from pathlib import Path
import pandas as pd
from . import BasicConfig, ScriptConfig, ResultDirConfig


@dataclass
class DataPreprocessingConfig:
    """Data preprocessing configuration class"""
    
    # Input file paths (must be provided during instantiation)
    profile: Union[str, Path, pd.DataFrame]
    phosfile: Union[str, Path, pd.DataFrame]
    sample_group: Union[str, Path, pd.DataFrame]
    mappingfile: Union[str, Path, pd.DataFrame]
    metadatafile: Union[str, Path, pd.DataFrame]

    # Optional input files
    prodiff: Optional[Union[str, Path, pd.DataFrame]] = None
    phosdiff: Optional[Union[str, Path, pd.DataFrame]] = None

    # Output directory
    outdir: Optional[Union[str, Path]] = None   # Path(BasicConfig.outdir +'/results/'+ ResultDirConfig.data_preprocessing)
    
    # Script path
    script_path: Optional[str] = None  # ScriptConfig.script_path
    
    # Omics names
    omics1_name: Optional[str] = None  # BasicConfig.omics1_name
    omics2_name: Optional[str] = None  # BasicConfig.omics2_name
    
    min_samples: ClassVar[int] = 3
    group_comparing: Optional[str] = None  # BasicConfig.group_comparing
    
    # Data filtering parameters (default values, can be overridden)
    pro_miss_value_ratio: float = 0.5
    phos_miss_value_ratio: float = 0.5

    # Missing value imputation parameters
    pro_imputation_method: Literal["zero","min","mean","median","knn","gmm","MLE"] = "knn"
    phos_imputation_method: Literal["zero","min","mean","median","knn","gmm","MLE"] = "knn"
    knn_k: ClassVar[int] = 5
    
    # Data normalization parameters
    pro_normalization_method: Literal["median","quantile","PQN","none"] = "median"
    phos_normalization_method: Literal["median","quantile","PQN","none"] = "median"
   
    def __post_init__(self):
        # outdir: use BasicConfig + ResultDirConfig if not provided
        if self.outdir is None:
            base = getattr(BasicConfig, "outdir", None)
            if base is None:
                # fallback to cwd if BasicConfig not set
                self.outdir = Path.cwd()
            else:
                self.outdir = Path(base) / getattr(ResultDirConfig, "data_preprocessing", "data_preprocessing")
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

        # group_comparing
        if not self.group_comparing:
            self.group_comparing = getattr(BasicConfig, "group_comparing", "T:N")

        # Enforce type consistency (for downstream usage)
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)

        # Ensure numeric default ranges are valid
        if not (0.0 <= float(self.pro_miss_value_ratio) <= 1.0):
            raise ValueError("pro_miss_value_ratio must be within [0,1]")
        if not (0.0 <= float(self.phos_miss_value_ratio) <= 1.0):
            raise ValueError("phos_miss_value_ratio must be within [0,1]")

    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return {
            # Input files
            'profile': str(cls.profile) if isinstance(cls.profile, (str, Path)) else 'DataFrame',
            'phosfile': str(cls.phosfile) if isinstance(cls.phosfile, (str, Path)) else 'DataFrame',
            'prodiff': str(cls.prodiff) if isinstance(cls.prodiff, (str, Path)) else 'DataFrame',
            'phosdiff': str(cls.phosdiff) if isinstance(cls.phosdiff, (str, Path)) else 'DataFrame',
            'sample_group': str(cls.sample_group) if isinstance(cls.sample_group, (str, Path)) else 'DataFrame',
            'metadatafile': str(cls.metadatafile) if isinstance(cls.metadatafile, (str, Path)) else 'DataFrame',
            
            # Output directory
            'outdir': str(cls.outdir),
            
            # Script path
            'script_path': cls.script_path,
            
            # Omics names
            'omics1_name': cls.omics1_name,
            'omics2_name': cls.omics2_name,
            
            # Identifier type
            'identified_type': cls.identified_type,
            
            # Data filtering parameters
            'pro_miss_value_ratio': cls.pro_miss_value_ratio,
            'phos_miss_value_ratio': cls.phos_miss_value_ratio,
            'pro_zero_value_ratio': cls.pro_zero_value_ratio,
            'phos_zero_value_ratio': cls.phos_zero_value_ratio,
            'pro_valid_value_ratio': cls.pro_valid_value_ratio,
            'phos_valid_value_ratio': cls.phos_valid_value_ratio,
            
            # Data normalization parameters
            'pro_normalization_method': cls.pro_normalization_method,
            'phos_normalization_method': cls.phos_normalization_method,
            'pro_scaling_method': cls.pro_scaling_method,
            'phos_scaling_method': cls.phos_scaling_method,
            
            # Data transformation parameters
            'pro_transformation_method': cls.pro_transformation_method,
            'phos_transformation_method': cls.phos_transformation_method,
            
            # Missing value imputation parameters
            'pro_imputation_method': cls.pro_imputation_method,
            'phos_imputation_method': cls.phos_imputation_method,
            'pro_imputation_value': cls.pro_imputation_value,
            'phos_imputation_value': cls.phos_imputation_value,
            'knn_k': cls.knn_k,
            
            # Batch correction parameters
            'batch_correction_method': cls.batch_correction_method,
            'batch_column': cls.batch_column,
            
            # Quality control parameters
            'qc_rsd_threshold': cls.qc_rsd_threshold,
            'qc_mad_threshold': cls.qc_mad_threshold,
            'remove_outsiders': cls.remove_outsiders,
            'outsider_sd_threshold': cls.outsider_sd_threshold
        }
    
    @classmethod
    def update_config(cls, **kwargs):
        """Update configuration parameters"""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
