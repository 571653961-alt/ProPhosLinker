#!/usr/bin/env python3
"""
ColorConfig: Comprehensive Color Palette for ProPhosLinker Visualizations

Centralized color management system for all visualization components in the 
proteomics/phosphoproteomics analysis pipeline. Organizes colors by analysis module 
with consistent two-tone scheme (orange: #ff7f00, blue: #012da7) and specialized 
palettes for networks, heatmaps, and functional enrichment.

Color Organization:
1. Pattern Analysis (PA, WGCNA) - Module identification colors
2. Differential Analysis - Quadrant, phosphorate, Mfuzz, networks
3. Functional Analysis - Enrichment gradients, KG networks, hub networks
4. Network Edge Types - 18 interaction types with semantic color coding

Usage:
    colors = ColorConfig.to_dict()           # All colors as dict
    groups = ColorConfig.get_color_groups()   # Grouped by analysis type
    ColorConfig.update_config(...)           # Dynamic updates
"""


from dataclasses import dataclass
from typing import ClassVar, List, Dict

@dataclass
class ColorConfig:
    """颜色配置类"""
    
    # # 1. Pattern analysis colors
    # #1.1 PA colors
    # PA_group1_color: ClassVar[str] = "#a03c32"
    # PA_group2_color: ClassVar[str] = "#1a5f6e"
    
    # # 1.3 WGCNA colors
    # WGCNA_pro_ME_color: ClassVar[str] = "#01344F"
    # WGCNA_phos_ME_color: ClassVar[str] = "#D12128"
    # WGCNA_pro_color: ClassVar[str] = "#a03c32"
    # WGCNA_phos_color: ClassVar[str] = "#43656C"
    # WGCNA_pheatmap_color: ClassVar[List[str]] = ["#43656C", "white", "#a03c32"]
    
    # 2. Differential analysis colors
    # 2.1 Omics differential quadrant colors
    # omics_diff_quadrant_pro_color: ClassVar[str] = "#a03c32"
    # omics_diff_quadrant_phos_color: ClassVar[str] = "#43656C"
    
    # # 2.2 Phosphorylation rate analysis colors
    # raw_pro_point_color: ClassVar[str] = '#01344F'
    # fitted_pro_point_color: ClassVar[str] = '#D12128'
    # fitted_pro_line_color: ClassVar[str] = '#01344F'
    # pro_fc_line_color: ClassVar[str] = '#4D613C'
    # total_phos_rate_line_color: ClassVar[str] = '#1a5f6e'
    # color_lst: ClassVar[List[str]] = ['#a03c32', '#1a5f6e', 'blue', 'orange', 'green', 'purple', 'cyan', 'yellow']
    
    # Mfuzz colors
    # Mfuzz_pluscolor: ClassVar[str] = "#a03c32"
    # Mfuzz_minuscolor: ClassVar[str] = '#1a5f6e'
    # Mfuzz_allcolor: ClassVar[str] = "grey60"
    # Mfuzz_onlyupcolor: ClassVar[str] = '#a03c32'
    # Mfuzz_onlydowncolor: ClassVar[str] = "#1a5f6e"
    # Mfuzz_updowncolor: ClassVar[str] = "#4D613C"
    # Mfuzz_notsigcolor: ClassVar[str] = "grey60"
    
    #   # 2.3 Differential network colors
    # differential_network_gradientn_color: ClassVar[List[str]] = ["#43656C", "white", "#a03c32"]
    # differential_network_edge_color_pos: ClassVar[str] = "#9b6a65"
    # differential_network_edge_color_neg: ClassVar[str] = "#5d8992"
    # differential_network_Enhanced_in_N: ClassVar[str] = "#5d8992"
    # differential_network_Enhanced_in_T: ClassVar[str] = "#9b6a65"
    # differential_network_Only_in_N: ClassVar[str] = "#0c2b32"
    # differential_network_Only_in_T: ClassVar[str] = "#381512"
    # differential_network_Conflict_relation : ClassVar[str] = "#808080"
    
    # # 3.1  Functional enrichment colors
    # omics_function_enrichment_color_gradient_low: ClassVar[str] = "#43656C"
    # omics_function_enrichment_color_gradient_high: ClassVar[str] = "#a03c32"
    
    # # 3.2 KG-based functional network node colors
    # kgbased_functional_network_node_up: ClassVar[str] = "#B83A2D"
    # kgbased_functional_network_node_down: ClassVar[str] = "#657f68"
    # kgbased_functional_network_node_nonsig: ClassVar[str] = "#E3C79F"
    # kgbased_functional_network_node_notdet: ClassVar[str] = "grey90"
    # kgbased_functional_network_node_disease_border_color: ClassVar[str] = "#535c54"
    # kgbased_functional_network_node_notdisease_border_color: ClassVar[str] = "#acabab"
    # kgbased_functional_network_function_enrich_color_gradient_low: ClassVar[str] = "#175663"
    # kgbased_functional_network_function_enrich_color_gradient_high: ClassVar[str] = "#90362d"
    
    # # 3.3 Hub protein network node colors
    # hub_protein_network_node_up: ClassVar[str] = "#B83A2D"
    # hub_protein_network_node_down: ClassVar[str] = "#657f68"
    # hub_protein_network_node_nonsig: ClassVar[str] = "#E3C79F"
    # hub_protein_network_node_notdet: ClassVar[str] = "grey90"
    # hub_protein_network_node_disease_border_color: ClassVar[str] = "#535c54"
    # hub_protein_network_node_notdisease_border_color: ClassVar[str] = "#acabab"
    # hub_protein_network_function_enrich_color_gradient_low: ClassVar[str] = "#175663"
    # hub_protein_network_function_enrich_color_gradient_high: ClassVar[str] = "#90362d"
    
    # Knowledge graph functional network edge type color mapping
    kgbased_functional_network_edge_type_colors: ClassVar[Dict[str, str]] = {
        "association": "#d7d6d6",           # Light gray
        "physical association": "#838181",  # Medium gray (same as binding)
        "binding": "#838181",               # Medium gray
        "direct interaction": "#d7d6d6",    # Light gray
        
        "activation": "#d62c0c",            # Orange-red
        "catalysis": "#7f00ff",             # Light purple
        "proximity": "#FF9301",             # Orange
        "reaction": "#114335",              # Dark green
        "phosphorylation": "#2FBE95",       # Teal
        "dephosphorylation": "#6B8E23",     # Olive green
        
        "ptmod": "#8C97D6",                 # Light purple
        "inhibition": "#0cb6d6",            # Cyan
        "expression": "#FCF402",            # Yellow
        "regulation": "#e4dadd",            # Purple-gray
        
        "colocalization": "#4c95cd",        # Sky blue
        "covalent binding": "#716F74",      # Violet
        "ubiquitination": "#FF4500",        # Orange-red
        "multiRel": "#9a6728"               # Brown
    }
    
    # Hub protein network edge type color mapping
    hub_protein_network_edge_type_colors: ClassVar[Dict[str, str]] = {
        "association": "#d7d6d6",           # Light gray
        "physical association": "#838181",  # Medium gray (same as binding)
        "binding": "#838181",               # Medium gray
        "direct interaction": "#d7d6d6",    # Light gray
        
        "activation": "#d62c0c",            # Orange-red
        "catalysis": "#7f00ff",             # Light purple
        "proximity": "#FF9301",             # Orange
        "reaction": "#114335",              # Dark green
        "phosphorylation": "#2FBE95",       # Teal
        "dephosphorylation": "#6B8E23",     # Olive green
        
        "ptmod": "#8C97D6",                 # Light purple
        "inhibition": "#0cb6d6",            # Cyan
        "expression": "#FCF402",            # Yellow
        "regulation": "#e4dadd",            # Purple-gray
        
        "colocalization": "#4c95cd",        # Sky blue
        "covalent binding": "#716F74",      # Violet
        "ubiquitination": "#FF4500",        # Orange-red
        "multiRel": "#9a6728"               # Brown
    }

    ############ test colors (active configuration)
    PA_group1_color: ClassVar[str] = "#ff7f00"
    PA_group2_color: ClassVar[str] = "#012da7"
    WGCNA_pro_ME_color: ClassVar[str] = "#FFCF14"
    WGCNA_phos_ME_color: ClassVar[str] = "#002FA7"
    WGCNA_pro_color: ClassVar[str] = "#ff7f00"
    WGCNA_phos_color: ClassVar[str] = "#012da7"
    WGCNA_pheatmap_color: ClassVar[List[str]] = ["#012da7", "white", "#ff7f00"]

    # 2. Differential analysis colors
    # 2.1 Omics differential quadrant colors
    omics_diff_quadrant_pro_color: ClassVar[str] = "#ff7f00"
    omics_diff_quadrant_phos_color: ClassVar[str] = "#012da7"
    # 2.2 Phosphorylation rate analysis colors
    raw_pro_point_color: ClassVar[str] = '#01344F'
    fitted_pro_point_color: ClassVar[str] = '#D12128'
    fitted_pro_line_color: ClassVar[str] = '#01344F'
    pro_fc_line_color: ClassVar[str] = '#4D613C'
    total_phos_rate_line_color: ClassVar[str] = "#ff7f00"
    color_lst: ClassVar[List[str]] = ["#ff7f00","#012da7",'blue','orange','green','purple','cyan','yellow']
    Mfuzz_minuscolor: ClassVar[str] = "#012da7"
    Mfuzz_pluscolor: ClassVar[str] = "#ff7f00"
    Mfuzz_allcolor: ClassVar[str] = "grey60"
    Mfuzz_onlyupcolor: ClassVar[str] = "#ff7f00"
    Mfuzz_onlydowncolor: ClassVar[str] = "#012da7"
    Mfuzz_updowncolor: ClassVar[str] = "#4D613C"
    Mfuzz_notsigcolor: ClassVar[str] = "grey60"
    # 2.3 Differential network colors
    differential_network_gradientn_color: ClassVar[List[str]] = ["#012da7", "white", "#ff7f00"]
    differential_network_edge_color_pos: ClassVar[str] = "#ffaa93"
    differential_network_edge_color_neg: ClassVar[str] = "#6090B8"
    differential_network_Enhanced_in_N: ClassVar[str] = "#41BBC8"
    differential_network_Enhanced_in_T: ClassVar[str] = "#77431E"
    differential_network_Only_in_N: ClassVar[str] = "#6090B8"
    differential_network_Only_in_T: ClassVar[str] = "#B7553C"
    differential_network_Conflict_relation: ClassVar[str] = "#808080"
    # 3.1 Functional enrichment colors
    omics_function_enrichment_color_gradient_low: ClassVar[str] = "#012da7"
    omics_function_enrichment_color_gradient_high: ClassVar[str] = "#ff7f00"
    # 3.2 KG-based functional network node colors
    kgbased_functional_network_node_up: ClassVar[str] = "#ff7f00"
    kgbased_functional_network_node_down: ClassVar[str] = "#012da7"
    kgbased_functional_network_node_nonsig: ClassVar[str] = "#E3C79F"
    kgbased_functional_network_node_notdet: ClassVar[str] = "grey90"
    kgbased_functional_network_node_disease_border_color: ClassVar[str] = "#B83A2D" 
    kgbased_functional_network_node_notdisease_border_color: ClassVar[str] = "#bec389"
    kgbased_functional_network_function_enrich_color_gradient_low: ClassVar[str] = "#012da7"
    kgbased_functional_network_function_enrich_color_gradient_high: ClassVar[str] = "#ff7f00"
    # 3.3 Hub protein network node colors
    hub_protein_network_node_up: ClassVar[str] = "#ff7f00"
    hub_protein_network_node_down: ClassVar[str] = "#012da7"
    hub_protein_network_node_nonsig: ClassVar[str] = "#E3C79F"
    hub_protein_network_node_notdet: ClassVar[str] = "grey90"
    hub_protein_network_node_disease_border_color: ClassVar[str] = "#B83A2D" 
    hub_protein_network_node_notdisease_border_color: ClassVar[str] = "#bec389"
    hub_protein_network_function_enrich_color_gradient_low: ClassVar[str] = "#012da7"
    hub_protein_network_function_enrich_color_gradient_high: ClassVar[str] = "#ff7f00"

    @classmethod
    def to_dict(cls) -> dict:
        """Convert to dictionary"""
        return {
            # Pattern analysis colors
            'PA_group1_color': cls.PA_group1_color,
            'PA_group2_color': cls.PA_group2_color,
            
            # WGCNA colors
            'WGCNA_pro_ME_color': cls.WGCNA_pro_ME_color,
            'WGCNA_phos_ME_color': cls.WGCNA_phos_ME_color,
            'WGCNA_pro_color': cls.WGCNA_pro_color,
            'WGCNA_phos_color': cls.WGCNA_phos_color,
            'WGCNA_pheatmap_color': cls.WGCNA_pheatmap_color,
            
            # Differential analysis colors
            'omics_diff_quadrant_pro_color': cls.omics_diff_quadrant_pro_color,
            'omics_diff_quadrant_phos_color': cls.omics_diff_quadrant_phos_color,
            
            # Phosphorylation rate analysis colors
            'raw_pro_point_color': cls.raw_pro_point_color,
            'fitted_pro_point_color': cls.fitted_pro_point_color,
            'fitted_pro_line_color': cls.fitted_pro_line_color,
            'pro_fc_line_color': cls.pro_fc_line_color,
            'total_phos_rate_line_color': cls.total_phos_rate_line_color,
            'color_lst': cls.color_lst,
            
            # Mfuzz colors
            'Mfuzz_pluscolor': cls.Mfuzz_pluscolor,
            'Mfuzz_minuscolor': cls.Mfuzz_minuscolor,
            'Mfuzz_allcolor': cls.Mfuzz_allcolor,
            'Mfuzz_onlyupcolor': cls.Mfuzz_onlyupcolor,
            'Mfuzz_onlydowncolor': cls.Mfuzz_onlydowncolor,
            'Mfuzz_updowncolor': cls.Mfuzz_updowncolor,
            'Mfuzz_notsigcolor': cls.Mfuzz_notsigcolor,
            
            # Differential network colors
            'differential_network_gradientn_color': cls.differential_network_gradientn_color,
            'differential_network_edge_color_pos': cls.differential_network_edge_color_pos,
            'differential_network_edge_color_neg': cls.differential_network_edge_color_neg,
            'differential_network_Enhanced_in_N': cls.differential_network_Enhanced_in_N,
            'differential_network_Enhanced_in_T': cls.differential_network_Enhanced_in_T,
            'differential_network_Only_in_N': cls.differential_network_Only_in_N,
            'differential_network_Only_in_T': cls.differential_network_Only_in_T,
            'differential_network_Conflict_relation': cls.differential_network_Conflict_relation,
            
            # Functional enrichment colors
            'omics_function_enrichment_color_gradient_low': cls.omics_function_enrichment_color_gradient_low,
            'omics_function_enrichment_color_gradient_high': cls.omics_function_enrichment_color_gradient_high,
            
            # KG-based functional network colors
            'kgbased_functional_network_node_up': cls.kgbased_functional_network_node_up,
            'kgbased_functional_network_node_down': cls.kgbased_functional_network_node_down,
            'kgbased_functional_network_node_nonsig': cls.kgbased_functional_network_node_nonsig,
            'kgbased_functional_network_node_notdet': cls.kgbased_functional_network_node_notdet,
            'kgbased_functional_network_node_disease_border_color': cls.kgbased_functional_network_node_disease_border_color,
            'kgbased_functional_network_node_notdisease_border_color': cls.kgbased_functional_network_node_notdisease_border_color,
            'kgbased_functional_network_function_enrich_color_gradient_low': cls.kgbased_functional_network_function_enrich_color_gradient_low,
            'kgbased_functional_network_function_enrich_color_gradient_high': cls.kgbased_functional_network_function_enrich_color_gradient_high,
            'kgbased_functional_network_edge_type_colors': cls.kgbased_functional_network_edge_type_colors,
            
            # Hub protein network colors
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
        """Update configuration parameters"""
        for key, value in kwargs.items():
            if hasattr(cls, key):
                setattr(cls, key, value)
            else:
                raise AttributeError(f"Invalid configuration parameter: {key}")
    
    @classmethod
    def get_color_groups(cls) -> Dict[str, List[str]]:
        """Get color configuration grouped by function"""
        return {
            'pattern_analysis': [
                'PA_group1_color', 'PA_group2_color',
                'WGCNA_pro_ME_color', 'WGCNA_phos_ME_color',
                'WGCNA_pro_color', 'WGCNA_phos_color',
                'WGCNA_pheatmap_color'
            ],
            'differential_analysis': [
                'omics_diff_quadrant_pro_color', 'omics_diff_quadrant_phos_color',
                'raw_pro_point_color', 'fitted_pro_point_color', 'fitted_pro_line_color',
                'pro_fc_line_color', 'total_phos_rate_line_color', 'color_lst',
                'Mfuzz_pluscolor', 'Mfuzz_minuscolor', 'Mfuzz_allcolor',
                'Mfuzz_onlyupcolor', 'Mfuzz_onlydowncolor', 'Mfuzz_updowncolor', 'Mfuzz_notsigcolor',
                'differential_network_gradientn_color', 'differential_network_edge_color_pos',
                'differential_network_edge_color_neg', 'differential_network_Enhanced_in_N',
                'differential_network_Enhanced_in_T', 'differential_network_Only_in_N',
                'differential_network_Only_in_T','differential_network_Conflict_relation'
            ],
            'functional_analysis': [
                'omics_function_enrichment_color_gradient_low', 'omics_function_enrichment_color_gradient_high',
                'kgbased_functional_network_node_up', 'kgbased_functional_network_node_down',
                'kgbased_functional_network_node_nonsig', 'kgbased_functional_network_node_notdet',
                'kgbased_functional_network_node_disease_border_color', 'kgbased_functional_network_node_notdisease_border_color',
                'kgbased_functional_network_function_enrich_color_gradient_low', 'kgbased_functional_network_function_enrich_color_gradient_high',
                'kgbased_functional_network_edge_type_colors',
                'hub_protein_network_node_up', 'hub_protein_network_node_down',
                'hub_protein_network_node_nonsig', 'hub_protein_network_node_notdet',
                'hub_protein_network_node_disease_border_color', 'hub_protein_network_node_notdisease_border_color',
                'hub_protein_network_function_enrich_color_gradient_low', 'hub_protein_network_function_enrich_color_gradient_high',
                'hub_protein_network_edge_type_colors'
            ]
        }