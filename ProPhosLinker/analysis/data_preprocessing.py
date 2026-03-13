#!/usr/bin/env python3
"""
ProPhosLinker: Data Preprocessing Module (0)

Comprehensive preprocessing pipeline for integrated proteomics and phosphoproteomics 
data analysis. This module handles the complete data cleaning and normalization 
workflow required for downstream differential analysis and network construction.

Key Features:
- Multi-format input support (TSV DataFrame)
- Sample identifier harmonization across omics datasets
- Flexible missing data filtering (feature-wise thresholds)
- 7 imputation methods: zero, min, mean, median, KNN, GMM, MLE (Gamma)
- 4 normalization methods: Median, Quantile, Z-score, PQN
- Outlier detection/removal and infinite value handling
- Automated validation and comprehensive reporting
- Group comparison validation (T:N format)
- Phosphosite-protein mapping validation

Workflow:
1. Input validation → Sample alignment → Missing data filtering
2. Outlier removal → Imputation → Normalization → Output generation
3. Comprehensive data quality reporting

"""

import re
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler
from scipy.stats import gamma
from sklearn.mixture import GaussianMixture

from ..config import DataPreprocessingConfig


class DataPreprocessing:
    """
    Comprehensive data preprocessing workflow for proteomics/phosphoproteomics integration.
    
    Responsibilities:
    - Validate input paths and parameter ranges
    - Load abundance matrices and metadata with format auto-detection
    - Harmonize sample identifiers across omics datasets
    - Filter features by missingness thresholds
    - Remove outliers and handle infinite values
    - Apply configurable imputation strategies
    - Perform normalization (median, quantile, z-score, PQN)
    - Generate quality reports and save downstream-ready outputs
    """
    
    
    def __init__(self, config: DataPreprocessingConfig):
        """
        Initialize preprocessing pipeline with configuration.
        
        Args:
            config: DataPreprocessingConfig with all input parameters
        """
        self.config = config
        
        def _load_input(x):
            if x is None:
                return None
            if isinstance(x, pd.DataFrame):
                return x.copy()
            p = Path(x)
            if not p.exists():
                raise FileNotFoundError(f"Input file not found: {p}")
            # assume TSV with index in first column; adjust if needed
            return pd.read_csv(p, sep='\t', index_col=0)
        # load required attributes so self.profile exists
        self.profile = _load_input(config.profile)
        self.phosfile = _load_input(config.phosfile)
        self.sample_group = _load_input(config.sample_group)
        self.mappingfile = _load_input(config.mappingfile)
        self.metadatafile = _load_input(config.metadatafile)
        self.prodiff = _load_input(config.prodiff)
        self.phosdiff = _load_input(config.phosdiff)
        # copy commonly used self.config values to instance attributes
        self.outdir = Path(config.outdir)
        self.script_path = config.script_path
        self.omics1_name = config.omics1_name
        self.omics2_name = config.omics2_name
        self.min_samples = config.min_samples
        self.group_comparing = config.group_comparing



    def validate_params(self):
        """Validate required file paths and high-level configuration values."""
        # 1) Required files must exist
        if not Path(self.config.profile).exists():
            raise FileNotFoundError(f"The total protein expression profile (-profile) file does not exist: {self.config.profile}")
        if not Path(self.config.phosfile).exists():
            raise FileNotFoundError(f"The phosphoprotein expression profile (-phosfile) file does not exist: {self.config.phosfile}")
        if not Path(self.config.mappingfile).exists():
            raise FileNotFoundError(f"The mappingfile profile (-mappingfile) file does not exist: {self.config.mappingfile}")
        if not Path(self.config.sample_group).exists():
            raise FileNotFoundError(f"The sample_group (-sample_group) file does not exist: {self.config.sample_group}")
        
        # 2) Validate group comparison format: "A:B"
        group_comparing_pattern = r'^[^:]+:[^:]+$'  
        if not bool(re.match(group_comparing_pattern, self.config.group_comparing)):
            raise ValueError("Parameter format error: must contain exactly one colon, with characters present both before and after it.")
        group_comparing_list = self.config.group_comparing.split(':')
        if len(group_comparing_list) != 2:
            raise ValueError("Parameter format error: must contain exactly one colon, with characters present both before and after it.")
        
        # 3) Validate numeric ranges
        if self.config.pro_miss_value_ratio< 0 or self.config.pro_miss_value_ratio> 1:
            raise ValueError("The pro_miss_value_ratiomust be within the range [0, 1].")
        if self.config.phos_miss_value_ratio< 0 or self.config.phos_miss_value_ratio> 1:
            raise ValueError("The phos_miss_value_ratiomust be within the range [0, 1].")
        
        # 4) Prepare output directory
        # Ensure outdir is a Path object
        if isinstance(self.config.outdir, str):
            self.config.outdir = Path(self.config.outdir)
        elif self.config.outdir is None:
            self.config.outdir = os.getcwd()
        elif not isinstance(self.config.outdir, Path):
            raise TypeError("outdir must be a string or Path object")
        if self.config.outdir.exists():
            if not os.access(self.config.outdir, os.W_OK):
                raise ValueError(f"The output directory already exists and overwriting is not allowed: {self.config.outdir}")
        # Create the output directory if it does not exist
        if not self.config.outdir.exists():
            self.config.outdir.mkdir(parents=True, exist_ok=True)
        try:
            self.config.outdir.mkdir(parents=True,exist_ok=True)
            print(f"✅ Output directory is ready: {self.config.outdir}", file=sys.stderr)
        except PermissionError:
            RuntimeError(f"No permission to create directory: {self.config.outdir}")
        except Exception as e:
            raise RuntimeError(f"Directory creation failed: {str(e)}")

    def check_load_data(self):
        """Load and validate all input datasets with format checking."""
        def validate_and_load_abundance(input_data):
            """Load abundance data with numeric conversion and format validation."""
            if isinstance(input_data, (str, Path)):
                sep = '\t' if str(input_data).endswith(('.tsv', '.txt')) else ','
                df = pd.read_csv(input_data, sep=sep, index_col=0, na_values=["NA", "-", "null", ""],low_memory=False)
            elif isinstance(input_data, pd.DataFrame):
                df = input_data.copy()
            else:
                raise TypeError("Input must be file path or DataFrame")
            
            # Validate index/column types
            if not all(isinstance(idx, str) for idx in df.index):
                raise ValueError("Row indices must be strings")
            if not all(isinstance(col, str) for col in df.columns):
                raise ValueError("Column names must be strings")
            return df.apply(pd.to_numeric, errors='coerce')
        
        def validate_and_load(input_data):
            if isinstance(input_data, (str, Path)):
                """Load metadata/annotation files."""
                sep = '\t' if str(input_data).endswith(('.tsv', '.txt')) else ','
                df = pd.read_csv(input_data, sep=sep, na_values=["NA", "-", "null", ""],low_memory=False)
            elif isinstance(input_data, pd.DataFrame):
                df = input_data.copy()
            else:
                raise TypeError("Error: Input must be a file path or a DataFrame")
            return df
        
        self.config.profile = validate_and_load_abundance(self.config.profile)
        self.config.phosfile = validate_and_load_abundance(self.config.phosfile)
        self.config.sample_group = validate_and_load(self.config.sample_group)
        if self.config.sample_group.columns[0] != 'Protein_sample':
            raise ValueError(f"Expected first column to be 'Protein_sample', but got '{self.config.sample_group.columns[0]}'")
        if self.config.sample_group.columns[1] != 'Phosphosite_sample':
            raise ValueError(f"Expected second column to be 'Phosphosite_sample', but got '{self.config.sample_group.columns[0]}'")
        if self.config.sample_group.columns[2] != 'Cor_sample':
            raise ValueError(f"Expected second column to be 'Cor_sample', but got '{self.config.sample_group.columns[0]}'")
        if self.config.sample_group.columns[3] != 'group':
            raise ValueError(f"Expected second column to be 'group', but got '{self.config.sample_group.columns[0]}'")
        sample_columns = ['Protein_sample', 'Phosphosite_sample', 'Cor_sample']
        for col in sample_columns:
            duplicated_mask = self.config.sample_group[col].duplicated()
            if duplicated_mask.any():
                duplicate_values = self.config.sample_group.loc[duplicated_mask, col].unique()
                raise ValueError(
                    f"Found duplicate values in column '{col}': {sorted(duplicate_values)}\n"
                    f"Total rows: {len(self.config.sample_group)}, Unique {col}: {self.config.sample_group[col].nunique()}"
                )

        self.config.mappingfile = validate_and_load(self.config.mappingfile)
        self.config.mappingfile.columns = [self.config.omics1_name,self.config.omics2_name]

        meta = self.config.metadatafile
        has_metadata = False
        if isinstance(meta, pd.DataFrame):
            has_metadata = not meta.empty
        elif isinstance(meta, (str, Path)):
            has_metadata = bool(str(meta).strip())
        elif meta is None:
            has_metadata = False
        else:
            has_metadata = False

        if has_metadata:
            self.config.metadatafile = validate_and_load(self.config.metadatafile)
            if 'sample' not in self.config.metadatafile.columns:
                raise ValueError(f"Required column 'sample' not found in metadata file. ")
            duplicate_samples = self.config.metadatafile['sample'].duplicated()
            if duplicate_samples.any():
                duplicate_values = self.config.metadatafile.loc[duplicate_samples, 'sample'].unique()
                duplicate_count = len(duplicate_values)
                raise ValueError(
                    f"Found {duplicate_count} duplicate sample(s) in metadatafile: {sorted(duplicate_values)}\n"
                    f"Total rows: {len(self.config.metadatafile)}, Unique samples: {self.config.metadatafile['sample'].nunique()}\n"
                    f"Duplicate entries: {self.config.metadatafile[self.config.metadatafile['sample'].isin(duplicate_values)].to_dict('records')}"
                )
        
        #############check_load_diff_file
        def check_load_diff_file(input_data,omics_name):
            if input_data is None:
                return None
            if isinstance(input_data, (str, Path)):
                sep = '\t' if str(input_data).endswith(('.tsv', '.txt')) else ','
                df = pd.read_csv(input_data, sep=sep, na_values=["NA", "-", "null", ""],low_memory=False)
            elif isinstance(input_data, pd.DataFrame):
                df = input_data.copy()
            else:
                raise TypeError("Input must be a file path (str) or a pandas DataFrame.")
            if df.shape[1] < 4:
                raise ValueError("Differential analysis file must contain at least 4 columns: ID, logFC, q_value, and class.")
            if df.isnull().any().any():
                raise ValueError("Differential analysis file must not contain any missing values.")
            df.columns = [omics_name,'logFC','adj.P.Val','class'] + list(df.columns[4:])
            return df

        self.config.prodiff = check_load_diff_file(self.config.prodiff,self.config.omics1_name)
        self.config.phosdiff = check_load_diff_file(self.config.phosdiff,self.config.omics2_name)
    
    def validate_group_sample(self):
        sample_group_unique = set(self.config.sample_group['group'].unique())
        group_comparing_set = set(self.config.group_comparing.split(':'))

        if sample_group_unique != group_comparing_set:
            raise ValueError(
                f"Group names mismatch! The sample_group file has groups: {sorted(sample_group_unique)}, "
                f"but comparing groups: {sorted(group_comparing_set)}"
            )

        sample_group_protein = set(self.config.sample_group['Protein_sample'])
        sample_group_phosphosite = set(self.config.sample_group['Phosphosite_sample'])
        sample_group_cor = set(self.config.sample_group['Cor_sample'])
        samples_profile = set(self.config.profile.columns)
        samples_phosfile = set(self.config.phosfile.columns)
        if sample_group_protein != samples_profile or sample_group_phosphosite != samples_phosfile:
            all_samples = sample_group_protein | samples_profile | samples_phosfile | sample_group_phosphosite
            missing_in_profile = sample_group_protein - samples_profile
            missing_in_phosfile = sample_group_phosphosite - samples_phosfile
            extra_in_profile = samples_profile - sample_group_protein
            extra_in_phosfile = samples_phosfile - sample_group_phosphosite
            
            error_msg = "Sample names mismatch across datasets!\n"
            if missing_in_profile:
                error_msg += f"Samples in group but missing in profile: {sorted(missing_in_profile)}\n"
            if missing_in_phosfile:
                error_msg += f"Samples in group but missing in phosfile: {sorted(missing_in_phosfile)}\n"
            if extra_in_profile:
                error_msg += f"Extra samples in profile: {sorted(extra_in_profile)}\n"
            if extra_in_phosfile:
                error_msg += f"Extra samples in phosfile: {sorted(extra_in_phosfile)}\n"
            
            raise ValueError(error_msg)

        protein_to_cor = dict(
            zip(
                self.config.sample_group['Protein_sample'],
                self.config.sample_group['Cor_sample']
            )
        )
        pro_new_columns = [protein_to_cor[col] for col in self.config.profile.columns]
        self.config.profile.columns = pro_new_columns
        Phosphosite_to_cor = dict(
            zip(
                self.config.sample_group['Phosphosite_sample'],
                self.config.sample_group['Cor_sample']
            )
        )
        phos_new_columns = [Phosphosite_to_cor[col] for col in self.config.phosfile.columns]
        self.config.phosfile.columns = phos_new_columns

        # Reorder the sample_group 'group' column to match the order in group_comparing (e.g., "T:N" keeps T before N)
        # 1. Determine the desired group order, e.g., ["T", "N"]
        group_order = self.config.group_comparing.split(':')

        # 2. Convert the 'group' column to an ordered categorical variable (defines the sort order)
        self.config.sample_group['group'] = pd.Categorical(
            self.config.sample_group['group'],
            categories=group_order,
            ordered=True
        )

        # 3. Key step: sort by the 'group' column so that T comes before N (preserving within-group order)
        self.config.sample_group = self.config.sample_group.sort_values('group').reset_index(drop=True)

        # 3. Check that all samples in metadatafile are present in sample_group
        # Determine metadatafile type to decide whether to load/validate
        meta = self.config.metadatafile
        has_metadata = False
        if isinstance(meta, pd.DataFrame):
            has_metadata = not meta.empty
        elif isinstance(meta, (str, Path)):
            has_metadata = bool(str(meta).strip())
        elif meta is None:
            has_metadata = False
        else:
            # Other types are considered invalid; skip or handle as needed
            has_metadata = False

        if has_metadata:
            metadata_samples = set(self.config.metadatafile['sample'])
            missing_samples = metadata_samples - sample_group_cor
            if missing_samples:
                raise ValueError(
                    f"Found {len(missing_samples)} sample(s) in metadatafile that are not present in sample_group: "
                    f"{sorted(missing_samples)}\n"
                    f"Metadata samples: {sorted(metadata_samples)}\n"
                    f"Sample group samples: {sorted(sample_group_cor)}")

        # Check that all indices in phosfile exist in the specified column of mappingfile
        phosfile_index = set(self.config.phosfile.index)
        mapping_omics = set(self.config.mappingfile[self.config.omics2_name])
        missing_indices = phosfile_index - mapping_omics
        if missing_indices:
            raise ValueError(
                f"Found {len(missing_indices)} index value(s) in profile that are not present in mappingfile['{self.config.omics2_name}']: "
                f"{sorted(missing_indices)[:10]}{'...' if len(missing_indices) > 10 else ''}\n"
                f"Profile indices: {len(phosfile_index)}, Mapping indices: {len(mapping_omics)}\n"
                f"Intersection: {len(phosfile_index & mapping_omics)}"
            )

    def filter_missing_data(self):
        # print(f"Shape before filtering missing values: {self.config.profile.shape}")
        self.config.profile = self.config.profile.loc[self.config.profile.isna().mean(axis=1) < self.config.pro_miss_value_ratio]
        self.config.phosfile = self.config.phosfile.loc[self.config.phosfile.isna().mean(axis=1) < self.config.phos_miss_value_ratio]
        # print(f"Shape after filtering missing values: {self.config.phosfile.shape}")

    def handle_outliers(self):
        def handle_outliers_basic(df):
            """Basic outlier handling: replace negative and infinite values with 0."""
            # Create a copy to avoid modifying the original data
            df_clean = df.copy()
            # Replace negative values with 0
            df_clean[df_clean < 0] = 0
            # Replace infinite values with 0
            df_clean = df_clean.replace([np.inf, -np.inf], 0)
            return df_clean
        self.config.profile = handle_outliers_basic(self.config.profile)
        self.config.phosfile = handle_outliers_basic(self.config.phosfile)

    def apply_imputation(self):
        """Apply missing data imputation."""
        def imputation(method, df):
            df_filled = df.copy()
            orig_mask = df.isna()  # Track original missing positions to ensure imputed values are non-negative

            # Implement different imputation methods
            if method == "zero":
                df_filled = df.fillna(0)

            elif method == "min":
                min_val = df.min().min()
                # Imputed values must be non-negative
                fill_value = max(min_val - 1 if pd.notnull(min_val) else 0, 0)
                df_filled = df.fillna(fill_value)

            elif method == "mean":
                # Fill each row (feature) with its mean; if the mean is negative, use 0
                def fill_mean(row):
                    mean_val = row.mean(skipna=True)
                    fill = max(mean_val if pd.notnull(mean_val) else 0, 0)
                    return row.fillna(fill)
                df_filled = df.apply(fill_mean, axis=1)

            elif method == "median":
                def fill_median(row):
                    med = row.median(skipna=True)
                    fill = max(med if pd.notnull(med) else 0, 0)
                    return row.fillna(fill)
                df_filled = df.apply(fill_median, axis=1)

            elif method == "knn":
                imputer = KNNImputer(n_neighbors=5, weights='distance')
                df_filled = pd.DataFrame(imputer.fit_transform(df), index=df.index, columns=df.columns)

            elif method == "gmm":
                df_filled = df.copy()
                for idx, row in df.iterrows():
                    mask = row.isnull()
                    if mask.any():
                        observed = row[~mask].values.reshape(-1, 1)
                        if len(observed) > 1:
                            gmm = GaussianMixture(n_components=2, random_state=0)
                            gmm.fit(observed)
                            fill_val = float(gmm.means_.min())
                            fill_val = max(fill_val, 0)  # ensure non-negative
                            row[mask] = fill_val
                            df_filled.loc[idx] = row

            elif method == "MLE":
                def gamma_impute(col):
                    data = pd.to_numeric(col.dropna(), errors='coerce')
                    data = data[np.isfinite(data)]
                    if data.empty or data.nunique() < 2:
                        fill_value = data.mean() if not data.empty else 0
                        fill_value = max(fill_value, 0)
                        col[col.isna()] = fill_value
                        return col
                    a, loc, scale = gamma.fit(data)
                    fill_size = col.isna().sum()
                    fill_values = gamma.rvs(a, loc, scale, size=fill_size)
                    fill_values = np.clip(fill_values, 0, None)  # ensure non-negative
                    col.loc[col.isna()] = fill_values
                    return col
                df_filled = df_filled.apply(gamma_impute)

                # Final step: clip imputed values at 0 only for the originally missing positions (some methods may produce negative values)
                if orig_mask.any().any():
                    # Clip only the cells that were originally NaN to avoid changing existing negative values
                    # (if you want to set all negatives to 0, you can use df_filled = df_filled.clip(lower=0))
                    df_filled = df_filled.copy()
                    # Clip the original missing positions column-by-column using a boolean mask
                    for col in df_filled.columns:
                        mask_col = orig_mask[col]
                        if mask_col.any():
                            df_filled.loc[mask_col, col] = df_filled.loc[mask_col, col].clip(lower=0)

            return df_filled

        if self.config.pro_imputation_method!= "none":
            self.config.profile = imputation(self.config.pro_imputation_method, self.config.profile)
        if self.config.phos_imputation_method!= "none":
            self.config.phosfile = imputation(self.config.phos_imputation_method, self.config.phosfile)

    def apply_normalization(self):
        def normalization(method,df):
            if method == "median":
                # Median centering: normalize by sample medians to adjust for sample differences
                df_qn = df.div(df.median(axis=0), axis=1)
            elif method == "quantile":
                # Quantile normalization
                rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
                df_qn = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
            elif method == "zscore":
                # Z-score normalization (standardization): normalize each column (sample) by subtracting the mean and dividing by the standard deviation (x - μ) / σ
                scaler = StandardScaler()
                df_qn = pd.DataFrame(
                    scaler.fit_transform(df),
                    index=df.index,
                    columns=df.columns
                )
            elif method == "PQN":
                # PQN (Probabilistic Quotient Normalization)
                # 1. Compute ratios of each sample (column) to a reference sample (typically the median sample)
                reference = df.median(axis=1)  # Use row (feature) medians as the reference profile
                quotients = df.div(reference, axis=0)  # Divide each sample by the reference profile
                factors = quotients.median(axis=0)    # Use the median ratio per sample as normalization factor
                df_qn = df.div(factors, axis=1)       # Normalize each sample by its own factor
            return df_qn

        if self.config.pro_normalization_method!= "none":
            self.config.profile = normalization(self.config.pro_normalization_method,self.config.profile)
        if self.config.phos_normalization_method!= "none":
            self.config.phosfile = normalization(self.config.phos_normalization_method,self.config.phosfile)

    def save_preprocessed_data(self):
        self.config.profile.index.name = self.config.omics1_name
        self.config.phosfile.index.name = self.config.omics2_name
        # Check self.config.profile and self.config.phosfile for missing and infinite values
        try:
            prof_arr = self.config.profile.to_numpy()
        except Exception:
            prof_arr = np.array(self.config.profile.values, dtype=float)

        try:
            phos_arr = self.config.phosfile.to_numpy()
        except Exception:
            phos_arr = np.array(self.config.phosfile.values, dtype=float)

        # counts
        profile_nan_total = int(self.config.profile.isna().sum().sum())
        phosfile_nan_total = int(self.config.phosfile.isna().sum().sum())

        profile_pos_inf = int(np.isposinf(prof_arr).sum())
        profile_neg_inf = int(np.isneginf(prof_arr).sum())
        phosfile_pos_inf = int(np.isposinf(phos_arr).sum())
        phosfile_neg_inf = int(np.isneginf(phos_arr).sum())

        # replace infinite values with NaN (cleaning)
        if profile_pos_inf + profile_neg_inf > 0:
            self.config.profile = self.config.profile.replace([np.inf, -np.inf], np.nan)
        if phosfile_pos_inf + phosfile_neg_inf > 0:
            self.config.phosfile = self.config.phosfile.replace([np.inf, -np.inf], np.nan)

        # per-column missing summary
        profile_col_missing = self.config.profile.isna().sum().rename('missing_count')
        phosfile_col_missing = self.config.phosfile.isna().sum().rename('missing_count')

        # save cleaned files
        self.config.profile.to_csv(self.config.outdir / Path(self.config.omics1_name.lower()+'_preprocessed.tsv'), sep='\t')
        self.config.phosfile.to_csv(self.config.outdir / Path(self.config.omics2_name.lower()+'_preprocessed.tsv'), sep='\t')

        # # also save per-column missing counts as separate tsvs for easy inspection
        # profile_col_missing.to_frame().to_csv(self.config.outdir / Path(self.config.omics1_name.lower()+'_missing_by_column.tsv'), sep='\t')
        # phosfile_col_missing.to_frame().to_csv(self.config.outdir / Path(self.config.omics2_name.lower()+'_missing_by_column.tsv'), sep='\t')

        # save sample_group / mapping / metadata as before
        self.config.sample_group = self.config.sample_group[["Cor_sample", "group"]]
        self.config.sample_group.columns = ['sample','group']
        self.config.sample_group.to_csv(self.config.outdir / Path('sample_list.tsv'), sep='\t', index=False)
        self.config.mappingfile.to_csv(self.config.outdir / Path('mapping_list.tsv'), sep='\t', index=False)
        # Decide whether to load/validate based on the metadatafile type
        meta = self.config.metadatafile
        has_metadata = False
        has_prodiff = False
        has_phosdiff = False
        def check_file_exists(file,file_name):
            has_file = False
            if file is None:
                has_file = False
            if isinstance(file, pd.DataFrame):
                has_file = not file.empty
            if isinstance(file, (str, Path)):
                has_file = bool(str(file).strip())
            if has_file:
                file.to_csv(self.config.outdir / Path(file_name), sep='\t', index=False)
        check_file_exists(self.config.metadatafile,'metadata.tsv')
        check_file_exists(self.config.prodiff,'prodiff.tsv')
        check_file_exists(self.config.phosdiff,'phosdiff.tsv')

        

        # final user message
        print(f'✅ Preprocessed files saved to: {self.config.outdir}', file=sys.stderr)

    def report_profile(self):
        """Expression profile overview."""
        pro_df = self.config.profile
        phos_df = self.config.phosfile
        omics1_name = self.config.omics1_name
        omics2_name = self.config.omics2_name
        mappingfile = self.config.mappingfile
        group_comparing = self.config.group_comparing
        sample_groups = self.config.sample_group
        meta_df = self.config.metadatafile
        outdir = self.config.outdir
        group1 = group_comparing.split(':')[0]
        group1_samples = sample_groups[sample_groups['group'] == group1]['sample']
        group2 = group_comparing.split(':')[1]
        group2_samples = sample_groups[sample_groups['group'] == group2]['sample']
        #samples info
        with open(outdir / Path('input_overview.txt'), 'w', encoding='utf-8') as file:
            file.write("---------------- Input Data Overview ----------------\n\n")
            file.write(f"\t{group1} VS. {group2}\n\n")
            file.write("Sample information:\n")
            file.write(f"\n\t{group1} group sample count: {len(group1_samples)}\n\n\t{(','.join(group1_samples))}\n")
            file.write(f"\n\t{group2} group sample count: {len(group2_samples)}\n\n\t{(','.join(group2_samples))}\n")
            file.write("\nProtein abundance table info:\n\n")
            detect_protein_count = pro_df.shape[0]
            detect_protein_name = pro_df.index.tolist()
        
            file.write(f"\tDetected {detect_protein_count} proteins\n\n")

            detect_phosphoproteinSite_count = phos_df.shape[0]
            # detect_phosphoprotein_name = mappingfile[mappingfile[omics2_name] == phos_df.index.tolist()][omics1_name]
            detect_phosphoprotein_name =  mappingfile[mappingfile[omics2_name].isin(phos_df.index)][omics1_name].unique()

            file.write(f"\tDetected {detect_phosphoproteinSite_count} phosphosites, {len(detect_phosphoprotein_name)} phosphoproteins, average {int(detect_phosphoproteinSite_count)/len(detect_phosphoprotein_name):.3f} phosphosites per phosphoprotein\n\n")
            # overlap protein and phosphoprotein
            overlapping_protein = [name for name in detect_protein_name if name in detect_phosphoprotein_name]
            file.write(f"\tDetected {len(overlapping_protein)} overlapping proteins between protein and phosphoprotein sets\n")
            if isinstance(meta_df, pd.DataFrame):
                has_metadata = not meta_df.empty
            elif isinstance(meta_df, (str, Path)):
                has_metadata = bool(str(meta_df).strip())
            elif meta_df is None:
                has_metadata = False
            else:
                # Other types are considered invalid; skip or handle as needed
                has_metadata = False

            if has_metadata:
                file.write("Phenotype information:\n")
                metadata_sample_name = meta_df['sample'].unique().tolist()
                file.write(f"\n\tNumber of samples with phenotype information: {len(metadata_sample_name)}\n\n\t{(','.join(metadata_sample_name))}\n")

        print(f'✅ See input overview: {outdir}\\input_overview.txt ', file=sys.stderr)




    def run(self):
        """
        Run the full preprocessing pipeline.

        Steps
        - Validate parameters and input paths
        - Load tables and perform schema checks
        - Harmonize samples and group labels
        - Filter by missingness, handle outliers, impute missing values
        - Normalize and write preprocessed outputs + a brief input overview
        """
        # 1) Parameter validation
        self.validate_params()
        # 2) Load and validate file contents
        self.check_load_data()
        # 3) Validate group/sample mapping and align sample IDs
        self.validate_group_sample()
        # 4) Filter features by missingness
        self.filter_missing_data()
        # 5) Handle outliers and non-finite values
        self.handle_outliers()
        # 6) Impute missing values
        self.apply_imputation()
        # 7) Normalize
        self.apply_normalization()
        # 8) Save preprocessed tables
        self.save_preprocessed_data()
        # 9) Write a brief input overview report
        self.report_profile()
