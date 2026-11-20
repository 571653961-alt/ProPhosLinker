#!/usr/bin/env python3
"""
Main analysis controller for the Protein-Phosphorylation Site Correlation Analysis (ProPhosLinker) pipeline.

This module exposes the high-level `main` function that orchestrates the core analysis workflow:
it interprets which steps to run, updates global configuration objects, builds per-step configuration
instances, checks Neo4j status, and then calls the concrete analysis classes.

Core responsibilities:
- Parse the `steps` parameter and decide which pipeline stages to execute
  (data_preprocessing, pattern_analysis, differential_analysis, functional_analysis).
- Convert user-specified fold-change thresholds from linear scale to log2 scale and
  propagate them into BasicConfig and related configuration objects.
- Construct configuration objects:
  * DataPreprocessingConfig
  * PatternAnalysisConfig
  * DifferentialConfig
  * FunctionalConfig
  using either defaults or an external configuration dictionary (`config_data`).

- Manage integration with Neo4j via Neo4jManager for knowledge-graph-based functional/network analysis.
- Call the analysis classes:
  * DataPreprocessing
  * PatternAnalysis
  * DifferentialAnalysis
  * FunctionalAnalysis
  in the correct order, depending on the selected steps.

Typical usage:
- This module is usually invoked from the CLI layer (e.g. `cli.py`), which parses command-line
  arguments and then forwards all resolved paths and parameters to `main(...)`.
"""

import sys
from pathlib import Path
import math

from .config import (
    BasicConfig, ScriptConfig, Neo4jConfig, ResultDirConfig,
    DataPreprocessingConfig, PatternAnalysisConfig, DifferentialConfig, FunctionalConfig
)
from .analysis import DataPreprocessing, PatternAnalysis, DifferentialAnalysis, FunctionalAnalysis
from .database import Neo4jManager


def parse_steps_param(steps_str):
    """Parse the steps parameter and return a dict of steps to run.
    
    Args:
        steps_str: Step string, either 'all' or a comma-separated list of step names,
                   e.g. 'data_preprocessing,differential_analysis'.
    
    Returns:
        dict: Boolean flags for each step, e.g. {'data_preprocessing': True, ...}.
    """
    valid_steps = {
        'data_preprocessing',
        'pattern_analysis',
        'differential_analysis',
        'functional_analysis'
    }
    
    # Initialize all steps to False
    steps_dict = {step: False for step in valid_steps}
    
    if steps_str.lower() == 'all':
        # Enable all steps
        steps_dict = {step: True for step in valid_steps}
    else:
        # Parse the specified steps
        requested_steps = [s.strip() for s in steps_str.split(',')]
        for step in requested_steps:
            if step in valid_steps:
                steps_dict[step] = True
            else:
                print(f"⚠️ Warning: Invalid step name: '{step}'")
                print(f"   Valid steps: {', '.join(sorted(valid_steps))}")
    
    return steps_dict


def main(pro_file, phos_file, sample_group, mapping_file,  
         metadata_file=None, group_comparing='T:N', outdir='./results',
         identified_type='SYMBOL', omics1_name="Protein", omics2_name="PhosProtein",
         pro_diff_file=None, phos_diff_file=None,
         pro_FC=2, pro_diff_q_val=0.05,
         phos_FC=2, phos_diff_q_val=0.05,
        disease=None,
         network_FC=2, neo4j_path=None, username='neo4j',
         password='neo4j', steps='all', config_data=None):
    
    # Parse steps parameter
    steps_dict = parse_steps_param(steps)
    print(f"\n📋 Analysis steps to be executed:")
    for step_name, is_enabled in steps_dict.items():
        status = "✅" if is_enabled else "⏭️"
        print(f"  {status} {step_name}")
    print()
    
    ########## Parameter update interface ###########
    pro_log2FC = math.log2(pro_FC)
    phos_log2FC = math.log2(phos_FC)
    network_log2FC = math.log2(network_FC)
    BasicConfig.update_config(disease=disease)
    BasicConfig.update_config(identified_type=identified_type)
    BasicConfig.update_config(outdir=Path(outdir).resolve())
    BasicConfig.update_config(group_comparing=group_comparing)
    BasicConfig.update_config(omics1_name=omics1_name)
    BasicConfig.update_config(omics2_name=omics2_name)
    BasicConfig.update_config(pro_log2FC=pro_log2FC)
    BasicConfig.update_config(phos_log2FC=phos_log2FC)
    BasicConfig.update_config(pro_diff_q_val=pro_diff_q_val)
    BasicConfig.update_config(phos_diff_q_val=phos_diff_q_val)
    BasicConfig.update_config(network_log2FC=network_log2FC)
    if disease:
        BasicConfig.update_config(disease=disease)
    if neo4j_path:
        Neo4jConfig.update_config(neo4j_path=neo4j_path)
    Neo4jConfig.update_config(username=username)
    Neo4jConfig.update_config(password=password)
    if pro_diff_file:
        prodiff_path = BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('prodiff.tsv')
    else:
        prodiff_path = None
    if phos_diff_file:
        phosdiff_path = BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('phosdiff.tsv')
    else:
        phosdiff_path = None
    ########## Parameter update interface ###########
    
    # Add project root directory to Python path
    sys.path.insert(0, Path(__file__).parent)
    # Create directory structure
    ResultDirConfig.create_directories(BasicConfig.outdir)

    if not config_data:
        # 0. Create data preprocessing config
        data_config = DataPreprocessingConfig(
            profile=pro_file,
            phosfile=phos_file,
            sample_group=sample_group,
            mappingfile=mapping_file,
            metadatafile=metadata_file,
            prodiff=pro_diff_file,
            phosdiff=phos_diff_file
        )
        if metadata_file:
            metadata_file = BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('metadata.tsv')
        # 1. Pattern analysis config
        pattern_config = PatternAnalysisConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            metadatafile=metadata_file,
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv')
        )
        # 2. Differential analysis config
        diff_config = DifferentialConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            mappingfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('mapping_list.tsv'),
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv'),
            prodiff=prodiff_path,
            phosdiff=phosdiff_path,
            differential_network_diff_pro_path=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics1_name + '.tsv'),
            differential_network_diff_phos_path=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv'),
            global_phosprotein_diff=str(
                BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis)
                / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv')
            ),
            plot_n_representative=5
        )
        # 3. Functional analysis config
        neo4j_mgr = Neo4jManager(
            username=Neo4jConfig.username,
            password=Neo4jConfig.password,
            neo4j_path=Neo4jConfig.neo4j_path
        )
        func_config = FunctionalConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            mappingfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('mapping_list.tsv'),
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv'),
            prodiff=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics1_name + '.tsv'),
            phosdiff=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv'),
            network_log2FC=BasicConfig.network_log2FC,
            username=Neo4jConfig.username,
            password=Neo4jConfig.password
        )
    else:
        data_config = DataPreprocessingConfig(
            profile=pro_file,
            phosfile=phos_file,
            sample_group=sample_group,
            mappingfile=mapping_file,
            metadatafile=metadata_file,
            prodiff=pro_diff_file,
            phosdiff=phos_diff_file,
            pro_miss_value_ratio=config_data['data_processing']['pro_miss_value_ratio'],
            phos_miss_value_ratio=config_data['data_processing']['phos_miss_value_ratio'],
            pro_imputation_method=config_data['data_processing']['pro_imputation_method'],
            phos_imputation_method=config_data['data_processing']['phos_imputation_method'],
            pro_normalization_method=config_data['data_processing']['pro_normalization_method'],
            phos_normalization_method=config_data['data_processing']['phos_normalization_method'],
        )
        if metadata_file:
            metadata_file = BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('metadata.tsv')
        # 1. Pattern analysis config
        pattern_config = PatternAnalysisConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv'),
            metadatafile=metadata_file,
            dim_rd_mtd=config_data['pa']['dim_rd_mtd'],
            NMF_pro_filter_num=config_data['nmf']['NMF_pro_filter_num'],
            NMF_phos_filter_num=config_data['nmf']['NMF_phos_filter_num'],
            WGCNA_pro_filter_num=config_data['wgcna']['WGCNA_pro_filter_num'],
            WGCNA_phos_filter_num=config_data['wgcna']['WGCNA_phos_filter_num'],
            protein_cor_method=config_data['wgcna']['protein_cor_method'],
            protein_corFun_tmp=config_data['wgcna']['protein_corFun_tmp'],
            protein_cluster_method=config_data['wgcna']['protein_cluster_method'],
            protein_networkType=config_data['wgcna']['protein_networkType'],
            protein_RsquareCut_val=config_data['wgcna']['protein_RsquareCut_val'],
            protein_mergingThresh=config_data['wgcna']['protein_mergingThresh'],
            protein_minModuleSize=config_data['wgcna']['protein_minModuleSize'],
            protein_SoftPower=config_data['wgcna']['protein_SoftPower'],
            phosphoprotein_cor_method=config_data['wgcna']['phosphoprotein_cor_method'],
            phosphoprotein_corFun_tmp=config_data['wgcna']['phosphoprotein_corFun_tmp'],
            phosphoprotein_cluster_method=config_data['wgcna']['phosphoprotein_cluster_method'],
            phosphoprotein_networkType=config_data['wgcna']['phosphoprotein_networkType'],
            phosphoprotein_RsquareCut_val=config_data['wgcna']['phosphoprotein_RsquareCut_val'],
            phosphoprotein_mergingThresh=config_data['wgcna']['phosphoprotein_mergingThresh'],
            phosphoprotein_minModuleSize=config_data['wgcna']['phosphoprotein_minModuleSize'],
            phosphoprotein_SoftPower=config_data['wgcna']['phosphoprotein_SoftPower'],
            module_cor_threshhold=config_data['wgcna']['module_cor_threshhold'],
            module_cor_p_adj=config_data['wgcna']['module_cor_p_adj'],
            module_cor_method=config_data['wgcna']['module_cor_method']
        )
        plot_style = (config_data.get('phosrate_subtyping', {}) or {}).get('plot', {}) or {}
        # 2. Differential analysis config
        diff_config = DifferentialConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            prodiff=prodiff_path,
            phosdiff=phosdiff_path,
            mappingfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('mapping_list.tsv'),
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv'),
            differential_network_diff_pro_path=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics1_name + '.tsv'),
            differential_network_diff_phos_path=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv'),
            robust_quantiles=config_data['phosrate_subtyping']['robust_quantiles'],
            high_overlap_coefficient=config_data['phosrate_subtyping']['high_overlap_coefficient'],
            high_overlap_iqr_ratio=config_data['phosrate_subtyping']['high_overlap_iqr_ratio'],
            partial_overlap_coefficient=config_data['phosrate_subtyping']['partial_overlap_coefficient'],
            fdr_threshold=config_data['phosrate_subtyping']['fdr_threshold'],
            effect_threshold_log2=config_data['phosrate_subtyping']['effect_threshold_log2'],
            conditional_spline_df=config_data['phosrate_subtyping']['conditional_spline_df'],
            prediction_grid_points=config_data['phosrate_subtyping']['prediction_grid_points'],
            support_tolerance=config_data['phosrate_subtyping']['support_tolerance'],
            # 2.2 conditional-curve figure style (config.yaml: phosrate_subtyping -> plot)
            plot_width=float(plot_style.get('width', 8.0)),
            plot_height=float(plot_style.get('height', 6.0)),
            plot_dpi=int(plot_style.get('dpi', 300)),
            plot_base_size=float(plot_style.get('base_size', 11.0)),
            plot_title_size=float(plot_style.get('title_size', 17.0)),
            plot_axis_title_size=float(plot_style.get('axis_title_size', 13.0)),
            plot_axis_text_size=float(plot_style.get('axis_text_size', 11.0)),
            plot_legend_size=float(plot_style.get('legend_size', 12.0)),
            plot_caption_size=float(plot_style.get('caption_size', 12.0)),
            plot_line_width=float(plot_style.get('line_width', 2.4)),
            plot_point_size=float(plot_style.get('point_size', 1.8)),
            plot_point_alpha=float(plot_style.get('point_alpha', 0.55)),
            plot_ribbon_alpha=float(plot_style.get('ribbon_alpha', 0.15)),
            plot_spine_width=float(plot_style.get('spine_width', 1.2)),
            plot_tick_width=float(plot_style.get('tick_width', 1.2)),
            plot_tick_length=float(plot_style.get('tick_length', 4.0)),
            plot_dash_color=str(plot_style.get('dash_color', 'grey45')),
            plot_dash_width=float(plot_style.get('dash_width', 1.1)),
            plot_panel_low=str(plot_style.get('panel_low', '#F7F9FC')),
            plot_panel_mid=str(plot_style.get('panel_mid', '#FFF9EF')),
            plot_panel_high=str(plot_style.get('panel_high', '#F9F5FA')),
            plot_caption_color=str(plot_style.get('caption_color', 'grey30')),
            plot_group_labels=str(plot_style.get('group_labels', 'N,T')),
            plot_n_representative=int(plot_style.get('n_representative', 5)),
            global_phosprotein_diff=str(
                BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis)
                / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv')
            ),
            differential_network_filter_num=config_data['differential_network']['differential_network_filter_num'],
            differential_network_FC_threshold=config_data['differential_network']['differential_network_FC_threshold'],
            differential_network_p_threshold=config_data['differential_network']['differential_network_p_threshold'],
            differential_network_p_value_type=config_data['differential_network']['differential_network_p_value_type'],
            differential_network_nBoots=config_data['differential_network']['differential_network_nBoots'],
            differential_network_bootnet_R_threshold=config_data['differential_network']['differential_network_bootnet_R_threshold'],
            differential_network_nCores=config_data['differential_network']['differential_network_nCores'],
            differential_network_stability_threshold=config_data['differential_network']['differential_network_stability_threshold'],
            differential_network_cor_method=config_data['differential_network']['differential_network_cor_method'],
            differential_network_edge_FC_threshold=config_data['differential_network']['differential_network_edge_FC_threshold'],
            differential_network_edge_p_threshold=config_data['differential_network']['differential_network_edge_p_threshold'],
            differential_network_max_subnet_num=config_data['differential_network']['differential_network_max_subnet_num'],
            differential_network_R_threshold=config_data['differential_network']['differential_network_R_threshold'],
            log_transform=config_data['differential_analysis']['log_transform'],
        )
        # 3. Functional analysis config
        neo4j_mgr = Neo4jManager(
            username=Neo4jConfig.username,
            password=Neo4jConfig.password,
            neo4j_path=Neo4jConfig.neo4j_path
        )
        func_config = FunctionalConfig(
            profile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics1_name.lower() + '_preprocessed.tsv'),
            phosfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path(BasicConfig.omics2_name.lower() + '_preprocessed.tsv'),
            mappingfile=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('mapping_list.tsv'),
            sample_group=BasicConfig.outdir / Path(ResultDirConfig.data_preprocessing) / Path('sample_list.tsv'),
            prodiff=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics1_name + '.tsv'),
            phosdiff=BasicConfig.outdir / Path(ResultDirConfig.omics_differential_analysis) / Path(BasicConfig.group_comparing.replace(':', '_vs_') + '_' + BasicConfig.omics2_name + '.tsv'),
            network_log2FC=BasicConfig.network_log2FC,
            GO_showCategory=config_data['functional_enrichment']['GO_showCategory'],
            KEGG_showCategory=config_data['functional_enrichment']['KEGG_showCategory'],
            kgbased_network_top_n=config_data['kgbased_network']['kgbased_network_top_n'],
            kgbased_network_top_nodes_visualization_num=config_data['kgbased_network']['kgbased_network_top_nodes_visualization_num'],
            kgbased_network_max_phosphoSite_displayed=config_data['kgbased_network']['kgbased_network_max_phosphoSite_displayed'],
            kgbased_network_node_filtering=config_data['kgbased_network']['kgbased_network_node_filtering'],
            kgbased_network_layout=config_data['kgbased_network']['kgbased_network_layout'],
            kgbased_network_comm_detection=config_data['kgbased_network']['kgbased_network_comm_detection'],
            hub_protein_top_n=config_data['hub_protein_network']['hub_protein_top_n'],
            hub_protein_d_num=config_data['hub_protein_network']['hub_protein_d_num'],
            hub_protein_network_layout=config_data['hub_protein_network']['hub_protein_network_layout'],
            username=Neo4jConfig.username,
            password=Neo4jConfig.password
        )

    ### Running pipeline ...
    if steps_dict['data_preprocessing']:
        print("🔄 Running: Data Preprocessing")
        DataPreprocessing(data_config).run()
    else:
        print("⏭️  Skipping: Data Preprocessing")
    
    if steps_dict['pattern_analysis']:
        print("\n🔄 Running: Pattern Analysis")
        PatternAnalysis(pattern_config).run()
    else:
        print("\n⏭️  Skipping: Pattern Analysis")
        
    if steps_dict['differential_analysis']:
        print("\n🔄 Running: Differential Analysis")
        DifferentialAnalysis(diff_config).run()
    else:
        print("\n⏭️  Skipping: Differential Analysis")
        
    if steps_dict['functional_analysis']:
        # Ensure database is running
        if neo4j_mgr.ensure_running():
            print("\n🔄 Running: Functional Analysis")
            print("✅ Neo4j is ready, starting work\n")
            FunctionalAnalysis(func_config).run()
        else:
            print("❌ Neo4j failed to start, please check manually!!!!!!!\n")
    else:
        print("\n⏭️  Skipping: Functional Analysis")

    print("\n✅ Analysis pipeline completed!")
