#!/usr/bin/env python3
"""
ScriptConfig: R Script Path Configuration for ProPhosLinker Pipeline

Centralized configuration for all R scripts used across the proteomics/phosphoproteomics 
analysis pipeline. Provides path management, validation utilities, and dictionary export 
for subprocess execution.

Script Organization:
1. Pattern Analysis (1 script)
2. Differential Analysis (3 scripts)  
3. Functional Analysis (5 scripts)
Total: 9 R scripts + base script path

Key Features:
- Inherits script_path from BasicConfig
- Full path construction (script_path + filename)
- Script existence validation
- Dictionary export for subprocess calls
- Key-based script lookup

Usage:
    paths = ScriptConfig.get_all_scripts()           # All script paths
    r_path = ScriptConfig.get_full_path('wgcna')     # Full R script path
    ScriptConfig.validate_paths()                    # Check all exist
"""


from dataclasses import dataclass
from typing import ClassVar, Union, Dict
from pathlib import Path
from .basic_config import BasicConfig


@dataclass
class ScriptConfig:
    """Script configuration class"""
    
    # Base script path
    script_path: ClassVar[str] = BasicConfig.script_path
    
    # Pattern analysis scripts
    pattern_analysis_rscript_path: ClassVar[Union[str, Path]] = '1.pattern_analysis.R'
    
    # Differential analysis scripts
    omics_differential_rscript_path: ClassVar[Union[str, Path]] = "2.1omics_differential_analysis.R"
    phosRate_quantile_subtyping_rscript_path: ClassVar[Union[str, Path]] = "2.2phosRate_quantile_subtyping.R"
    differential_network_rscript_path: ClassVar[Union[str, Path]] = "2.3differential_network.R"
    
    # Functional analysis scripts
    omics_function_enrichment_rscript_path: ClassVar[Union[str, Path]] = '3.1functional_enrichment.R'
    kgbased_functional_network_rscript_path: ClassVar[Union[str, Path]] = '3.2functionnal_interaction_network_visualization.R'
    kgbased_functional_network_community_rscript_path: ClassVar[Union[str, Path]] = '3.2functionnal_interaction_network_community_visualisation.R'
    hub_protein_functional_network_rscript_path: ClassVar[Union[str, Path]] = "3.3network_hubgen_visualization.R"
    function_enrich_rscript_path: ClassVar[Union[str, Path]] = '3.2functional_enrichment_function.R'
    
    @classmethod
    def get_all_scripts(cls) -> Dict[str, str]:
        """Get all script configurations"""
        return {
            'script_path': cls.script_path,
            'pattern_analysis_rscript_path': cls.pattern_analysis_rscript_path,
            'omics_differential_rscript_path': cls.omics_differential_rscript_path,
            'phosRate_quantile_subtyping_rscript_path': cls.phosRate_quantile_subtyping_rscript_path,
            'differential_network_rscript_path': cls.differential_network_rscript_path,
            'omics_function_enrichment_rscript_path': cls.omics_function_enrichment_rscript_path,
            'kgbased_functional_network_rscript_path': cls.kgbased_functional_network_rscript_path,
            'kgbased_functional_network_community_rscript_path': cls.kgbased_functional_network_community_rscript_path,
            'hub_protein_functional_network_rscript_path': cls.hub_protein_functional_network_rscript_path,
            'function_enrich_rscript_path': cls.function_enrich_rscript_path
        }
    
    @classmethod
    def get_full_path(cls, script_key: str) -> Path:
        """Get full script path"""
        scripts = cls.get_all_scripts()
        if script_key not in scripts:
            raise KeyError(f"Script key '{script_key}' not found")
        
        script_name = scripts[script_key]
        return Path(cls.script_path) / script_name
    
    @classmethod
    def validate_paths(cls) -> bool:
        """Validate script paths exist"""
        import os
        missing_scripts = []
        
        for key, script_name in cls.get_all_scripts().items():
            if key == 'script_path':
                continue
                
            full_path = cls.get_full_path(key)
            if not os.path.exists(full_path):
                missing_scripts.append((key, str(full_path)))
        
        if missing_scripts:
            print("Missing scripts:")
            for key, path in missing_scripts:
                print(f"  {key}: {path}")
            return False
        return True
    
    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return cls.get_all_scripts()
