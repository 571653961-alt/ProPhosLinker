#!/usr/bin/env python3
"""
BasicConfig: Core Configuration Class for ProPhosLinker Pipeline

Centralized configuration management for the entire proteomics/phosphoproteomics 
analysis pipeline. Defines default paths, thresholds, omics naming, and analysis 
parameters using dataclass with ClassVar for type safety and immutability.

Key Features:
- Script paths (auto-detected relative to module location)
- Output directory management
- Omics naming conventions (Protein/PhosProtein)
- Differential analysis thresholds (log2FC, q-values)
- Identifier type (UNIPROT/SYMBOL)
- Disease association settings
- Network analysis parameters

Usage:
    config = BasicConfig.to_dict()  # Get as dictionary
    BasicConfig.update_config(outdir='./my_results')  # Update parameters
"""

from dataclasses import dataclass
from typing import ClassVar, Union, Literal
from pathlib import Path


@dataclass
class BasicConfig:
    """Basic configuration class"""
    
    # Base script path
    script_path: ClassVar[str] = Path(__file__).resolve().parent.parent / "scripts"
    
    # Output directory
    outdir: ClassVar[Union[str, Path]] = './results'

    # Report resource path
    resource_dir: ClassVar[Union[str, Path]] = Path(__file__).resolve().parent.parent / 'resources'
    
    # Omics names
    omics1_name: ClassVar[str] = 'Protein'
    omics2_name: ClassVar[str] = 'PhosProtein'
    
    # Group comparison
    group_comparing: ClassVar[str] = 'Tumor:Normal'
    
    # Differential analysis thresholds
    pro_log2FC: ClassVar[float] = 1
    phos_log2FC: ClassVar[float] = 1
    pro_diff_q_val: ClassVar[float] = 0.05
    phos_diff_q_val: ClassVar[float] = 0.05
    phosRate_log2FC: ClassVar[float] = 1
    
    # Identifier type
    identified_type: ClassVar[Literal['UNIPROT', 'SYMBOL']] = 'UNIPROT'
    
    # Disease
    disease: ClassVar[str] = ''  # 'pancreatic ductal adenocarcinoma'

    # Network differential analysis
    network_log2FC: ClassVar[float] = 1
    
    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return {
            'script_path': cls.script_path,
            'outdir': str(cls.outdir),
            'omics1_name': cls.omics1_name,
            'omics2_name': cls.omics2_name,
            'group_comparing': cls.group_comparing,
            'pro_log2FC': cls.pro_log2FC,
            'phos_log2FC': cls.phos_log2FC,
            'pro_diff_q_val': cls.pro_diff_q_val,
            'phos_diff_q_val': cls.phos_diff_q_val,
            'phosRate_log2FC': cls.phosRate_log2FC,
            'identified_type': cls.identified_type,
            'disease': cls.disease,
            'network_log2FC': cls.network_log2FC
        }
    
    @classmethod
    def update_config(cls, **kwargs):
        """Update configuration parameters"""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
