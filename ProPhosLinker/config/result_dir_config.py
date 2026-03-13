#!/usr/bin/env python3
"""
ResultDirConfig: Standardized Result Directory Structure

Dataclass defining the complete hierarchical directory structure for ProPhosLinker 
analysis pipeline outputs. Provides consistent organization across all analysis 
modules with main categories and specialized subdirectories for reproducibility.

Directory Hierarchy:
0. Data_Preprocessing/           <- Raw → Preprocessed data
1. Pattern_Analysis/             <- Procrustes, Subtyping, WGCNA
2. Differential_Analysis/        <- Omics DE, PhosRate, Stable Networks
3. Functional_Analysis/          <- Enrichment, KG Networks, Hub Networks

Key Features:
- 4 main analysis categories
- 9 specialized subdirectories  
- Automatic directory creation
- Path construction utilities
- Dictionary export for configuration

Usage:
    ResultDirConfig.create_directories('./results')  # Create all dirs
    path = ResultDirConfig.get_full_path('./results', 'wgcna')
"""


from dataclasses import dataclass
from typing import ClassVar, Dict
from pathlib import Path


@dataclass
class ResultDirConfig:
    """Result directory configuration class"""
    
    # Base directory structure
    data_preprocessing: ClassVar[str] = "0.Data_Preprocessing"
    pattern_analysis: ClassVar[str] = "1.Pattern_Analysis"
    differential_analysis: ClassVar[str] = "2.Differential_Analysis"
    functional_analysis: ClassVar[str] = "3.Functional_Analysis"
    
    # Subdirectories
    procrustes: ClassVar[str] = "1.Pattern_Analysis/1.1Procrustes"
    molecular_subtyping: ClassVar[str] = "1.Pattern_Analysis/1.2Molecular_Subtyping"
    wgcna: ClassVar[str] = "1.Pattern_Analysis/1.3WGCNA"
    omics_differential_analysis: ClassVar[str] = "2.Differential_Analysis/2.1Omics_Differential_Analysis"
    phosrate_subtyping: ClassVar[str] = "2.Differential_Analysis/2.2PhosRate_Subtyping"
    stable_differential_network: ClassVar[str] = "2.Differential_Analysis/2.3Stable_Differential_Network"
    omics_function_enrichment: ClassVar[str] = "3.Functional_Analysis/3.1Omics_Function_Enrichment"
    kgbased_functional_network: ClassVar[str] = "3.Functional_Analysis/3.2KGbased_Functional_Network"
    hub_protein_functional_network: ClassVar[str] = "3.Functional_Analysis/3.3Hub_Protein_Functional_Network"
    
    @classmethod
    def get_all_dirs(cls) -> Dict[str, str]:
        """Get all directory configurations"""
        return {
            'data_preprocessing': cls.data_preprocessing,
            'pattern_analysis': cls.pattern_analysis,
            'differential_analysis': cls.differential_analysis,
            'functional_analysis': cls.functional_analysis,
            'procrustes': cls.procrustes,
            'molecular_subtyping': cls.molecular_subtyping,
            'wgcna': cls.wgcna,
            'omics_differential_analysis': cls.omics_differential_analysis,
            'phosrate_subtyping': cls.phosrate_subtyping,
            'stable_differential_network': cls.stable_differential_network,
            'omics_function_enrichment': cls.omics_function_enrichment,
            'kgbased_functional_network': cls.kgbased_functional_network,
            'hub_protein_functional_network': cls.hub_protein_functional_network
        }
    
    @classmethod
    def get_full_path(cls, base_dir: str, *sub_dirs: str) -> Path:
        """Get full path"""
        path = Path(base_dir)
        for sub_dir in sub_dirs:
            path = path / sub_dir
        return path
    
    @classmethod
    def create_directories(cls, base_output_dir: str):
        """Create all result directories"""
        base_path = Path(base_output_dir)
        
        for dir_name in cls.get_all_dirs().values():
            dir_path = base_path / dir_name
            dir_path.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {dir_path}")
    
    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return cls.get_all_dirs()
