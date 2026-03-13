#!/usr/bin/env python3
"""
Command-line interface for the Protein-Phosphorylation Site Correlation Analysis (ProPhosLinker) pipeline.

This script defines the CLI entry point `pro-phos-cor`, parses user inputs with argparse,
validates required files, optionally loads a YAML configuration file, prints a summary of
all effective parameters, and then calls the core analysis function in `ProPhosLinker.main`.

Main argument groups:
- Required parameters:
  * --pro_file / -p        : Protein expression data file (CSV/TSV format)
  * --phos_file / -s       : Phosphoprotein site data file (CSV/TSV format)
  * --sample_group / -g    : Sample grouping file
  * --mapping_file / -m    : Protein–phosphorylation site mapping file

- Common parameters:
  * --pro_diff_file        : Protein differential analysis result file
  * --phos_diff_file       : Phosphoprotein differential analysis result file
  * --neo4j_path           : Neo4j installation path
  * --outdir / -o          : Output directory (default: ./results)
  * --metadata_file        : Sample/clinical metadata file
  * --group_comparing / -c : Comparison groups, e.g. "T:N"
  * --identified_type / -i : Identifier type (UNIPROT or SYMBOL)
  * --pro_FC, --phos_FC    : Fold-change thresholds for protein and phosphosite
  * --pro_diff_q_val       : q-value threshold for protein differential analysis
  * --phos_diff_q_val      : q-value threshold for phosphoprotein differential analysis
  * --phosRate_FC          : Phosphorylation rate fold-change threshold
  * --disease              : Disease type / project description
  * --network_FC           : log2FC threshold for network analysis
  * --steps                : Analysis steps to run (all, data_preprocessing, pattern_analysis,
                             differential_analysis, functional_analysis)

- System configuration:
  * --config               : YAML configuration file
  * --username             : Neo4j username (default: neo4j)
  * --password             : Neo4j password (default: neo4j)

Usage examples are provided in the argparse epilog.
"""

import argparse
import sys
import os
import yaml
import subprocess


def main():
    """Command-line entry function."""
    parser = argparse.ArgumentParser(
        description="Protein-Phosphorylation Site Correlation Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Basic usage
            pro-phos-cor \
                --pro_file "protein_abundance.tsv" \
                --phos_file "phosphoprotein_abundance.tsv" \
                --sample_group "compare_groups.tsv" \
                --mapping_file "protein_phosphoproSite.tsv" \
                --group_comparing 'T:N' \
                --outdir "." \
                --config "config.yaml" \
                --pro_FC 2 \
                --phos_FC 2 \
                --network_FC 2 \
                --password "neo4j_password"

            # Use configuration file
                pro-phos-cor \
                    --pro_file "protein_abundance.tsv" \
                    --phos_file "phosphoprotein_abundance.tsv" \
                    --sample_group "compare_groups.tsv" \
                    --mapping_file "protein_phosphoproSite.tsv" \
                    --metadata_file "clinical_table_140.tsv" \
                    --group_comparing 'T:N' \
                    --outdir "." \
                    --config "config.yaml" \
                    --identified_type 'SYMBOL' \
                    --disease "pancreatic ductal adenocarcinoma" \
                    --pro_FC 2 \
                    --phos_FC 2 \
                    --network_FC 2 \

                    --password "neo4j_password"
                    """
    )
    # ==========================================================================
    # 1. Basic required parameters
    # ==========================================================================
    required_group = parser.add_argument_group('Required parameters')
    required_group.add_argument(
        '--pro_file', '-p', required=True, type=str,
        help='Protein data file (CSV/TSV format)'
    )
    required_group.add_argument(
        '--phos_file', '-s', required=True, type=str,
        help='Phosphoprotein site data file (CSV/TSV format)'
    )
    required_group.add_argument(
        '--sample_group', '-g', required=True, type=str,
        help='Samples and group data file (CSV/TSV format)'
    )
    required_group.add_argument(
        '--mapping_file', '-m', required=True, type=str,
        help='Protein–phosphorylation site mapping data file (CSV/TSV format)'
    )

    # ==========================================================================
    # 2. Common optional parameters
    # ==========================================================================
    basic_group = parser.add_argument_group('Common parameters')
    basic_group.add_argument(
        '--pro_diff_file', type=str, default=None,
        help='Protein differential analysis result file'
    )
    basic_group.add_argument(
        '--phos_diff_file', type=str, default=None,
        help='Phosphoprotein differential analysis result file'
    )
    basic_group.add_argument(
        '--neo4j_path', type=str, default=None,
        help='Neo4j installation path'
    )
    basic_group.add_argument(
        '--outdir', '-o', type=str, default='./results',
        help='Output directory (default: ./results)'
    )
    basic_group.add_argument(
        '--metadata_file', type=str, default=None,
        help='Metadata file (CSV/TSV format)'
    )
    basic_group.add_argument(
        '--group_comparing', '-c', type=str, default='T:N',
        help='Comparison groups (default: T:N)'
    )
    basic_group.add_argument(
        '--identified_type', '-i', type=str, default='SYMBOL',
        choices=['UNIPROT', 'SYMBOL'],
        help='Identifier type (default: SYMBOL)'
    )
    basic_group.add_argument(
        '--pro_FC', type=float, default=2.0,
        help='Protein FC threshold (default: 2.0)'
    )
    basic_group.add_argument(
        '--phos_FC', type=float, default=2.0,
        help='Phosphoprotein FC threshold (default: 2.0)'
    )
    basic_group.add_argument(
        '--omics1_name', type=str, default="Protein",
        help='Name for omics 1 (default: Protein)'
    )
    basic_group.add_argument(
        '--omics2_name', type=str, default="PhosProtein",
        help='Name for omics 2 (default: PhosProtein)'
    )
    basic_group.add_argument(
        '--pro_diff_q_val', type=float, default=0.05,
        help='Protein differential analysis q-value threshold (default: 0.05)'
    )
    basic_group.add_argument(
        '--phos_diff_q_val', type=float, default=0.05,
        help='Phosphoprotein differential analysis q-value threshold (default: 0.05)'
    )
    basic_group.add_argument(
        '--phosRate_FC', type=float, default=2.0,
        help='Phosphorylation rate FC threshold (default: 2.0)'
    )
    basic_group.add_argument(
        '--disease', type=str, default=None,
        help='Disease type / project description'
    )
    basic_group.add_argument(
        '--network_FC', type=float, default=2.0,
        help='Network analysis log2FC threshold (default: 2.0)'
    )
    basic_group.add_argument(
        '--steps', type=str, default='all',
        help=(
            'Steps to run (default: all). Options: all, data_preprocessing, '
            'pattern_analysis, differential_analysis, functional_analysis. '
            'Separate multiple steps with comma (e.g., "pattern_analysis,differential_analysis")'
        )
    )

    # ==========================================================================
    # 3. System configuration parameters
    # ==========================================================================
    system_group = parser.add_argument_group('System configuration')
    system_group.add_argument(
        '--config', type=str, default=None,
        help='Configuration file path (YAML format)'
    )
    system_group.add_argument(
        '--username', type=str, default='neo4j',
        help='Neo4j username (default: neo4j)'
    )
    system_group.add_argument(
        '--password', type=str, default='neo4j',
        help='Neo4j password (default: neo4j)'
    )

    # ==========================================================================
    # 4. Other options
    # ==========================================================================
    parser.add_argument(
        '--version', '-v', action='version',
        version='ProPhosLinker 0.1.0'
    )

    # Parse arguments
    args = parser.parse_args()

    # Validate required parameters
    if not validate_required_args(args):
        sys.exit(1)

    # Load configuration file (optional, ignore failure)
    config_data = load_config_file(args.config) if args.config else {}

    # Print argument summary
    print_args_summary(args, config_data)

    # Metadata file is optional; nothing to do if not provided
    if not args.metadata_file:
        args.metadata_file

    # Auto-detect neo4j bin if not provided (Windows `where neo4j`)
    if not args.neo4j_path:
        args.neo4j_path = (
            os.path.dirname(
                subprocess.run(
                    'where neo4j',
                    shell=True,
                    capture_output=True,
                    text=True
                ).stdout.strip().split('\n')[0]
            )
            if subprocess.run('where neo4j', shell=True, capture_output=True).returncode == 0
            else None
        )
        print(f"Neo4j bin path: {args.neo4j_path}")

    # Call main analysis function
    try:
        from ProPhosLinker.main import main as analysis_main

        # Pass all parameters to main analysis function
        analysis_main(
            pro_file=args.pro_file,
            phos_file=args.phos_file,
            sample_group=args.sample_group,
            mapping_file=args.mapping_file,
            metadata_file=args.metadata_file,
            group_comparing=args.group_comparing,
            outdir=args.outdir,
            identified_type=args.identified_type,
            omics1_name=args.omics1_name,
            omics2_name=args.omics2_name,
            pro_FC=args.pro_FC,
            phos_FC=args.phos_FC,
            pro_diff_q_val=args.pro_diff_q_val,
            phos_diff_q_val=args.phos_diff_q_val,
            phosRate_FC=args.phosRate_FC,
            disease=args.disease,
            network_FC=args.network_FC,
            neo4j_path=args.neo4j_path,
            pro_diff_file=args.pro_diff_file,
            phos_diff_file=args.phos_diff_file,
            username=args.username,
            password=args.password,
            steps=args.steps,
            config_data=config_data  # pass configuration data
        )

    except ImportError as e:
        print(f"❌ Failed to import analysis module: {e}")
        print("Please ensure the ProPhosLinker package is correctly installed")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Analysis execution failed: {e}")
        sys.exit(1)


def validate_required_args(args):
    """Validate required input files and Neo4j path (if provided)."""
    required_files = [
        (args.pro_file, 'Protein expression data file'),
        (args.phos_file, 'Phosphorylation data file'),
        (args.sample_group, 'Sample grouping file'),
        (args.mapping_file, 'Protein–phosphorylation site mapping file')
    ]

    all_files_exist = True
    for file_path, desc in required_files:
        if not os.path.exists(file_path):
            print(f"❌ {desc} does not exist: {file_path}")
            all_files_exist = False

    # Only check neo4j_path if provided; if not, warn but do not abort
    if args.neo4j_path:
        if not os.path.exists(args.neo4j_path):
            print(f"❌ Neo4j path does not exist: {args.neo4j_path}")
            all_files_exist = False
    else:
        print(
            "⚠️ Neo4j path not specified. Some features (such as database import/query) "
            "may not work. You can specify the path with --neo4j_path, or detect and set it later."
        )

    return all_files_exist


def load_config_file(config_path):
    """Load configuration file and return a dict; return empty dict on failure."""
    if not config_path or not os.path.exists(config_path):
        print(f"⚠️ Configuration file does not exist: {config_path}, using default configuration")
        return {}

    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)

        if config_data is None:
            print(f"⚠️ Configuration file is empty: {config_path}, using default configuration")
            return {}

        print(f"✅ Configuration file loaded successfully: {config_path}")
        return config_data

    except yaml.YAMLError as e:
        print(f"⚠️ Configuration file format error: {e}, using default configuration")
        return {}
    except Exception as e:
        print(f"⚠️ Failed to read configuration file: {e}, using default configuration")
        return {}


def print_args_summary(args, config_data):
    """Print a summary of CLI arguments and configuration sections."""
    print("=" * 70)
    print("📊 ProPhosLinker analysis argument summary")

    print("📁 Input files:")
    print(f"  Protein data: {args.pro_file}")
    print(f"  Phosphorylation data: {args.phos_file}")
    print(f"  Sample grouping: {args.sample_group}")
    print(f"  Mapping file: {args.mapping_file}")
    if args.metadata_file:
        print(f"  Metadata file: {args.metadata_file}")

    print("\n⚙️ Analysis configuration:")
    print(f"  Output directory: {args.outdir}")
    print(f"  Comparison groups: {args.group_comparing}")
    print(f"  Identifier type: {args.identified_type}")
    print(
        f"  Differential thresholds: pro_FC={args.pro_FC}, phos_FC={args.phos_FC}, "
        f"pro_diff_q_val={args.pro_diff_q_val}, phos_diff_q_val={args.phos_diff_q_val}"
    )

    print("\n🔧 System configuration:")
    print(f"  Neo4j path: {args.neo4j_path}")

    if config_data:
        config_sections = list(config_data.keys())
        print(f"\n📋 Configuration sections: {', '.join(config_sections)}")
    else:
        print("\n📋 Using default configuration parameters")

    print("=" * 70)
    print("🚀 Start analysis...")


if __name__ == '__main__':
    main()
