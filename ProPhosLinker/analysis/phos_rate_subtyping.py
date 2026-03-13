#!/usr/bin/env python3
"""
ProPhosLinker: Phosphorylation-rate quantile subtyping module for the ProPhosLinker pipeline.

This module defines the `PhosRateSubtyping` class and `phos_rate_subtyping()` function 
(part of ProPhosLinker toolkit), implementing a novel quantile-based phosphorylation-rate 
analysis workflow for subtype discovery.

Core workflow (4-step pipeline):
1) **Quantile transformation & fitting**: For each protein, sort samples by protein intensity 
   quantiles and fit a smooth protein intensity curve using linear interpolation (outlier removal).

2) **Phosphorylation rate computation**: Calculate phosphorylation-rate log2FC as:
PhosRate_log2FC(q) = [PhosSite(q) - Pro(q)]_G1 - [PhosSite(q) - Pro(q)]_G2

text
where q = quantile rank across group samples.

3) **Signal denoising**: Apply LOESS smoothing (statsmodels) to remove noise from quantile profiles.

4) **Mfuzz clustering**: Call Rscript for fuzzy clustering of denoised phosphorylation-rate profiles 
across quantile positions, identifying conserved patterns.

Key innovations:
- **Quantile normalization across groups** avoids batch effects in direct group comparisons
- **Outlier-aware curve fitting** with IQR-based outlier detection and interpolation
- **Comprehensive validation** of inputs, parameters, sample consistency
- **Traditional vs. quantile comparison** to identify sites missed by classical methods
- **Rich visualization** with fitted curves, quantile-rate profiles, and top/bottom diff sites

Typical usage:
- Used by `DifferentialAnalysis.phosrate_subtyping()` (step 2.2) after omics differential analysis.
- Generates `fitted_phosRate_quantile.tsv`, denoised profiles, Mfuzz clusters, and diagnostic plots.
"""


from typing import Union
from pathlib import Path
import pandas as pd
import numpy as np
import os
import sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.interpolate import interp1d
import subprocess
import gc

class PhosRateSubtyping:
    def __init__(
        self,
        *,
        profile : Union[str, Path, pd.DataFrame],
        phosphoprofile : Union[str, Path, pd.DataFrame],
        phosphopro_pro_cor : Union[str, Path, pd.DataFrame],
        prodiff : Union[str, Path, pd.DataFrame],
        phosphoprodiff : Union[str, Path, pd.DataFrame],
        #omics name
        omics1_name : str = 'Pro',
        omics2_name : str = 'Phos',
        group_vs : str = 'T:N',
        group_info:Union[str, Path, pd.DataFrame],
        outdir : Union[str, Path],
        r_script_path :Union[str, Path] = "E:/pro_phosphpro/script/phosRate_quantile_subtyping.R",
        phos_rate_FC: float = 1.0,
        quantile_rank : int = 100,
        quantile_rank_step : int = 1,
        frac_param: float = 1,
        cluster_num: int = 9,
        plot_sites_num : int = 7,
        plot_site_ids : list = None,
        vs_total_phos_rate_upper: float = 0.7,
        vs_total_phos_rate_lower: float = 0.4,
        top_tail_diff_site_num : int = 5,
        #
        raw_pro_point_color : str='#01344F',
        fitted_pro_point_color : str= '#D12128',
        fitted_pro_line_color: str = '#01344F',
        pro_fc_line_color : str = '#4D613C',
        total_phos_rate_line_color : str = '#1a5f6e',
        color_lst : list = ['#a03c32','#1a5f6e','blue','orange','green','purple','cyan','yellow'],
        Mfuzz_pluscolor : str = "#a03c32",
        Mfuzz_minuscolor : str = '#1a5f6e',
        Mfuzz_allcolor : str = "grey60",
        Mfuzz_onlyupcolor : str = '#a03c32',
        Mfuzz_onlydowncolor : str = "#1a5f6e",
        Mfuzz_updowncolor : str = "#4D613C",
        Mfuzz_notsigcolor : str = "grey60",
        

        ):
        #input file：
        self.profile = profile
        self.phosphoprofile = phosphoprofile
        self.phosphopro_pro_cor = phosphopro_pro_cor
        self.prodiff = prodiff
        self.phosphoprodiff = phosphoprodiff
        #
        self.omics1_name = omics1_name
        self.omics2_name = omics2_name
        self.group_vs = group_vs
        self.group_info = group_info

        #output directory
        self.outdir = outdir

        #analytical parameter
        self.phos_rate_FC = phos_rate_FC
        self.quantile_rank = quantile_rank
        self.quantile_rank_step = quantile_rank_step
        self.vs_total_phos_rate_upper = vs_total_phos_rate_upper
        self.vs_total_phos_rate_lower = vs_total_phos_rate_lower

        #plot parameter
        self.plot_sites_num = plot_sites_num
        self.plot_site_ids = plot_site_ids
        self.top_tail_diff_site_num = top_tail_diff_site_num
        # colors
        self.raw_pro_point_color = raw_pro_point_color
        self.fitted_pro_point_color = fitted_pro_point_color
        self.fitted_pro_line_color = fitted_pro_line_color
        self.color_lst = color_lst
        self.pro_fc_line_color = pro_fc_line_color
        self.total_phos_rate_line_color = total_phos_rate_line_color
        #
        self.Mfuzz_minuscolor = Mfuzz_minuscolor
        self.Mfuzz_pluscolor = Mfuzz_pluscolor
        self.Mfuzz_allcolor = Mfuzz_allcolor
        self.Mfuzz_onlyupcolor = Mfuzz_onlyupcolor
        self.Mfuzz_onlydowncolor = Mfuzz_onlydowncolor
        self.Mfuzz_updowncolor = Mfuzz_updowncolor
        self.Mfuzz_notsigcolor = Mfuzz_notsigcolor
        #denoise parameter
        self.frac_param = frac_param
        #cluster parameter
        self.cluster_num = cluster_num

        self.r_script_path = r_script_path


        self._validate_params()
        self._data_load()


    def _validate_params(self):
        # 1.Ensure outdir is a Path object
        if isinstance(self.outdir, str):
            self.outdir = Path(self.outdir)
        elif self.outdir is None:
            self.outdir = os.getcwd()
        elif not isinstance(self.outdir, Path):
            raise TypeError("outdir must be a string or Path object")
        if self.outdir.exists():
            if self.outdir.is_file():
                raise ValueError(f"Output path is a file, not a directory: {self.outdir}")
            if not os.access(self.outdir, os.W_OK):
                raise ValueError(f"Output directory exists and is not writable: {self.outdir}")
        
        # Create the output directory if it does not exist
        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)
        try:
            self.outdir.mkdir(parents=True,exist_ok=True)
            print(f"✅ Output directory ready: {self.outdir}", file=sys.stderr)
        except PermissionError:
            raise RuntimeError(f"No permission to create directory: {self.outdir}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory: {str(e)}")
    
        # 2.Ensure all input files are Path objects or DataFrames
        for attr in ['profile', 'phosphoprofile', 'phosphopro_pro_cor', 'group_info']:
            value = getattr(self, attr)
            if isinstance(value, str):
                setattr(self, attr, Path(value))
            elif not isinstance(value, (Path, pd.DataFrame)):
                raise TypeError(f"{attr} must be a string, Path object, or DataFrame")
            # Ensure the input files exist if they are Path objects
            if isinstance(value, Path) and not value.exists():
                raise FileNotFoundError(f"{value} does not exist")
            # Ensure the input files are readable if they are Path objects
            if isinstance(value, Path) and not value.is_file():
                raise ValueError(f"{value} is not a valid file")
            # Ensure the input files are DataFrames if they are not Path objects
            if isinstance(value, pd.DataFrame):
                if value.empty:
                    raise ValueError(f"{attr} DataFrame is empty")
                if not all(col in value.columns for col in ['ID', 'Name']):
                    raise ValueError(f"{attr} DataFrame must contain 'ID' and 'Name' columns")

        # 3.parameter check
        if self.phos_rate_FC<= 0 :
            raise ValueError("phos rate FC (log2FC) must be > 0")
        if self.quantile_rank<= 0 :
            raise ValueError("quantile_rank must be > 0")
        if self.quantile_rank_step<= 0 :
            raise ValueError("quantile_rank_step must be > 0")
        if self.quantile_rank <= self.quantile_rank_step:
            raise ValueError("quantile_rank must be greater than quantile_rank_step")
        if self.quantile_rank_step <= 0:
            raise ValueError("quantile_rank_step must be > 0")
        if self.plot_sites_num <=0 or self.plot_sites_num >7:
            raise ValueError("plot_sites_num must be an integer in (0, 7]")
        if self.frac_param <=0 or self.frac_param >1:
            raise ValueError("denoise frac_param must be a float in (0, 1]")
        if self.cluster_num <=1 or self.cluster_num >20:
            raise ValueError("cluster_num must be an integer in (1, 20]")
        if self.vs_total_phos_rate_upper <= self.vs_total_phos_rate_lower:
            raise ValueError("vs_total_phos_rate_upper must be greater than vs_total_phos_rate_lower")
        if self.vs_total_phos_rate_lower <0 or self.vs_total_phos_rate_lower >1:
            raise ValueError("vs_total_phos_rate_lower must be in [0, 1]")
        if self.vs_total_phos_rate_upper <0 or self.vs_total_phos_rate_upper >1:
            raise ValueError("vs_total_phos_rate_upper must be in [0, 1]")
        if not isinstance(self.group_vs, str) or ':' not in self.group_vs:
            raise ValueError("group_vs format error: expected 'group1:group2', e.g. 'T:N'")
        if self.top_tail_diff_site_num <=0 or self.top_tail_diff_site_num >2000:
            raise ValueError("top_tail_diff_site_num must be an integer in (0, 2000]")
    def _data_load(self):
        def datafile_load_validation(file,type):
            if not isinstance(file, pd.DataFrame):
                    # Load data from a tab-delimited file
                    if(type == 'mumber'):
                        df = pd.read_csv(file, sep='\t', index_col=0, na_values=["NA", "-", "null", ""])
                        df = df.apply(pd.to_numeric, errors='coerce')
                    else:
                        df = pd.read_csv(file,sep='\t',header=0)
            else:
                df = file
            if df.empty:
                raise ValueError(file.as_posix()+" expression profile data is empty")
            if df.index.duplicated().any():
                raise ValueError(file.as_posix()+" has duplicated row indices")
            return df
        
         # log2+1 transformation (safely handle common placeholders and avoid invalid values)
        def safe_log2p1(df):
            df_num = df.replace(['_Inf', '-Inf', 'Inf', 'inf', '-inf', 'NA', 'nan'], np.nan).apply(pd.to_numeric, errors='coerce')
            # If values < -1 exist, clip to -0.999999 (or set to NaN if preferred)
            if (df_num < -1).any().any():
                print(f"Warning: found values < -1 in data; clipping to -0.999999 before log2p1", file=sys.stderr)
                df_num = df_num.clip(lower=-0.999999)
            return np.log2(df_num + 1)
        
        self.profile = datafile_load_validation(self.profile,'mumber')
        # Apply log2+1 transformation to profile if needed (skip if data is already log2/normalized)
        try:
            self.profile = safe_log2p1(self.profile)
        except Exception as e:
            raise RuntimeError(f"log2+1 transform failed for profile: {e}")

        self.phosphoprofile = datafile_load_validation(self.phosphoprofile,'mumber')
        # Apply log2+1 transformation to phosphoprofile if needed (skip if data is already log2/normalized)
        try:
            self.phosphoprofile = safe_log2p1(self.phosphoprofile)
        except Exception as e:
            raise RuntimeError(f"log2+1 transform failed for phosphoprofile: {e}")

        self.phosphopro_pro_cor = datafile_load_validation(self.phosphopro_pro_cor,'string')
        self.group_info = datafile_load_validation(self.group_info,'string')
        self.prodiff = datafile_load_validation(self.prodiff,'mumber')
        self.phosphoprodiff = datafile_load_validation(self.phosphoprodiff,'mumber')

        min_group_sample_num = min(len(self.group_info[self.group_info['group'] == self.group_info['group'].unique()[0]]),len(self.group_info[self.group_info['group'] == self.group_info['group'].unique()[1]]))
        if self.quantile_rank <min_group_sample_num:
            raise ValueError("To avoid overfitting, quantile_rank must be greater than the smallest group sample size: " + str(min_group_sample_num))

        for group_name_vs in self.group_vs.split(':'):
            if group_name_vs not in self.group_info['group'].unique().tolist():
                raise ValueError(f"Group '{group_name_vs}' from group_vs not found in group_info")
        for group_name in self.group_info['group'].unique().tolist():
            if group_name not in self.group_vs.split(':'):
                print(f"⚠️ Warning: group '{group_name}' is not included in group_vs '{self.group_vs}'", file=sys.stderr)
        # Ensure all samples in group_info are present in profile and phosphoprofile
        profile_samples = set(self.profile.columns)
        phosphoprofile_samples = set(self.phosphoprofile.columns)
        group_info_samples = set(self.group_info['sample'])
        missing_in_profile = group_info_samples - profile_samples
        missing_in_phosphoprofile = group_info_samples - phosphoprofile_samples
        if missing_in_profile:
            raise ValueError(f"The following samples from group_info are missing in profile: {', '.join(missing_in_profile)}")
        if missing_in_phosphoprofile:
            raise ValueError(f"The following samples from group_info are missing in phosphoprofile: {', '.join(missing_in_phosphoprofile)}")



def pro_group_trans_intensity(config,quantile_rank,group_name,protein_name):
    # Target protein intensity
    pro_intensity = config.profile.loc[[protein_name]].T
    # print(pro_intensity)
    # Target protein phosphorylation site intensities
    phos_sites = config.phosphopro_pro_cor[config.phosphopro_pro_cor[config.omics1_name] == protein_name][config.omics2_name]
    phos_nomissing_sites = [site for site in phos_sites if site in config.phosphoprofile.index]
    phos_intensity = config.phosphoprofile.loc[phos_nomissing_sites].T
    # print(phos_intensity.head())
    # Combine protein and phospho-site intensities
    combined_pro_phos_df = pd.concat([pro_intensity,phos_intensity],axis=1)
    # print(combined_pro_phos_df.head())

    # Group-based fitting
    # For group 1
    group_sample_names = config.group_info[config.group_info['group'] == group_name]['sample']
    # print(group_sample_names)
    pro_intenstity_group = combined_pro_phos_df.loc[group_sample_names]
    # print(pro_intenstity_group.head())
    # add rank quantile_rank
    pro_intenstity_sorted_group= pro_intenstity_group.sort_values(by=protein_name)
    # print(pro_intenstity_sorted_group)
    pro_intenstity_sorted_group['rank'] = list(range(1,len(pro_intenstity_sorted_group)+1))
    pro_intenstity_sorted_group['quantile_rank'] = pro_intenstity_sorted_group['rank'] / len(pro_intenstity_sorted_group) *100
    # print(pro_intenstity_sorted_group.head())

    # Handle outliers in quantile_rank vs pro_intensity
    x_data = pro_intenstity_sorted_group['quantile_rank']
    y_data = pro_intenstity_sorted_group[protein_name]
    Q1 = y_data.quantile(0.25)
    Q3 = y_data.quantile(0.75)
    IQR = Q3 - Q1

    # Define outlier boundaries
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    # Create boolean mask to identify non-outliers
    non_outlier_mask = (y_data >= lower_bound) & (y_data <= upper_bound)

    x_clean = x_data[non_outlier_mask]
    y_clean = y_data[non_outlier_mask]

    # Interpolation fitting
    pro_intensity_interp = interp1d(
                x_clean,
                y_clean,
                kind='linear',
                fill_value='extrapolate'
            )
    
    pro_fitted = pro_intensity_interp(quantile_rank)

    trans_df = pd.DataFrame({
        'quantile_rank':quantile_rank,
        'intensity_fitted':pro_fitted
    })
    
    site_intensity_dict = {}
    for site_id in phos_nomissing_sites:
        target_quantiles = quantile_rank
        site_near_intensity = []
        for q in target_quantiles:
            distances = np.abs(pro_intenstity_sorted_group['quantile_rank'] - q)
            nearest_idx = distances.idxmin()
            nearest_site_intensity = pro_intenstity_sorted_group.loc[nearest_idx][site_id]
            site_near_intensity.append(nearest_site_intensity)
        site_intensity_dict[site_id] = site_near_intensity

    # Merge at once
    site_intensity_df = pd.DataFrame(site_intensity_dict, index=trans_df.index)
    trans_df = pd.concat([trans_df, site_intensity_df], axis=1)
        
    return trans_df,pro_intenstity_sorted_group

# group fitted result plot
def fitted_result_plot(config,pro_intenstity_sorted_group,trans_df,plot_site_ids,protein_name,group_name,ax):
    ax.scatter(pro_intenstity_sorted_group['quantile_rank'],
                pro_intenstity_sorted_group[protein_name],
                color=config.raw_pro_point_color,
                label='Raw pro intensity')
    ax.scatter(
            trans_df['quantile_rank'],
            trans_df['intensity_fitted'],
            color=config.fitted_pro_point_color,
            alpha=0.4,
            label='Fitted pro intensity'
        )
    ax.plot(
            trans_df['quantile_rank'],
            trans_df['intensity_fitted'],
            color=config.fitted_pro_line_color,
            linewidth=2,
            label=' Fitted pro intensity trend'
            )
    #phos sites
    sites_id = plot_site_ids[:config.plot_sites_num]
    for i in range(0,len(sites_id)):
        ax.scatter(pro_intenstity_sorted_group['quantile_rank'],
                    pro_intenstity_sorted_group[sites_id[i]],
                    color=config.color_lst[i],
                    alpha=0.2,
                    # label=sites_id[i].replace(protein_id+'_',protein_name+':')+' intensity'),
                    label = sites_id[i]+"'s intensity"
                    )
        
    for i in range(0,len(sites_id)):
        ax.scatter(trans_df['quantile_rank'],
                    trans_df[sites_id[i]],
                    color=config.color_lst[i],
                    # label=sites_id[i].replace(protein_id+'_',protein_name+':')+' nearest intensity '
                    label=sites_id[i]+"'s nearest intensity"
                    )
    
    ax.set_ylabel(protein_name +" and Phosphorylation Sites Intensity ("+group_name+")",fontsize=16)
    
# quantile_rank - phos site plot
def quantile_phos_site_plot(config,adjusted_pro_site,protein_name,plot_site_ids,total_phos_rate,ax):
    ax.plot(
            adjusted_pro_site.index,
            adjusted_pro_site[config.omics1_name+'_log2FC'],
            color= config.pro_fc_line_color,
            linewidth = 5,
            alpha=0.6,
            label=protein_name +"'s log2FC")
    sites_id = plot_site_ids[:config.plot_sites_num]
    for i in range(0,len(sites_id)):
        site = sites_id[i]
        ax.plot(
            adjusted_pro_site.index,
            adjusted_pro_site[site],
            color=config.color_lst[i],
            linewidth = 2,
            alpha=0.9,
            label=site+' '+config.omics1_name+' rate log2FC')
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.axhline(y=config.phos_rate_FC, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.axhline(y=-config.phos_rate_FC, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    if total_phos_rate:
        ax.axhline(y=total_phos_rate, color=config.total_phos_rate_line_color, linestyle='--', linewidth=2, alpha=0.8,label='total '+config.omics2_name+' rate log2FC: '+str(round(total_phos_rate,3)))

    ax.set_xlabel('Quantile',fontsize=16)
    ax.set_ylabel('log2FC',fontsize=16)
    # ax.set_title(protein_name +"'s pro and phos rate log2FC ("+config.group_vs.replace(':',' vs. ')+")")

def protein_phos_rate_fitted_analysis(config,protein_name,plot_tage=False,total_phos_rate=False):
    group_names = config.group_vs.split(':')
    quantile_rank = np.arange(0,config.quantile_rank + 1,config.quantile_rank_step)
    plot_site_ids = config.plot_site_ids
    group1_pro_trans_intensity,group1_raw_pro_intenstity_sorted = pro_group_trans_intensity(config,quantile_rank,group_names[0],protein_name)
    group2_pro_trans_intensity,group2_raw_pro_intenstity_sorted = pro_group_trans_intensity(config,quantile_rank,group_names[1],protein_name)

    adjusted_pro_site = pd.DataFrame({
        config.omics1_name+'_log2FC':(group1_pro_trans_intensity['intensity_fitted']-group2_pro_trans_intensity['intensity_fitted']),
    },index=group1_pro_trans_intensity.index)
    # Collect all new columns into a dict
    new_cols = {}
    for site in group1_pro_trans_intensity.columns.tolist()[2:]:
        new_cols[site] = group1_pro_trans_intensity[site] - group2_pro_trans_intensity[site] - adjusted_pro_site[config.omics1_name+'_log2FC']

    # Merge at once
    new_cols_df = pd.DataFrame(new_cols, index=adjusted_pro_site.index)
    adjusted_pro_site = pd.concat([adjusted_pro_site, new_cols_df], axis=1)

    if(not plot_site_ids):
        plot_site_ids = group1_pro_trans_intensity.columns.tolist()[2:]
    
    if plot_tage:
        fig, axes = plt.subplots(3, 1, figsize=(18, 18))  # 3 rows x 1 column

            # First plot
        fitted_result_plot(config, group1_raw_pro_intenstity_sorted, group1_pro_trans_intensity, plot_site_ids, protein_name, group_names[0], axes[0])
        axes[0].legend(loc='lower left')

        # Second plot
        fitted_result_plot(config, group2_raw_pro_intenstity_sorted, group2_pro_trans_intensity, plot_site_ids, protein_name, group_names[1], axes[1])
        axes[1].legend(loc='lower left')

        # Third plot
        quantile_phos_site_plot(config, adjusted_pro_site, protein_name, plot_site_ids, total_phos_rate, axes[2])
        # axes[2].legend()

        for ax in axes:
            ax.legend(
                loc='center left',
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0,
                fontsize=14  # Set label font size
            )

        fig.suptitle(protein_name +"'s Quantile - Phosphorylation Rate Analysis ("+config.group_vs.replace(':',' vs. ')+")", fontsize=22, y=0.99)
        plt.tight_layout()
        plt.savefig(str(config.outdir)+"/"+protein_name +"-"+plot_site_ids[0]+"_Quantile_Phosphorylation_Rate.png")
        plt.close(fig)
        gc.collect()
        
    return adjusted_pro_site

def get_phosRate_quantile_df(config):
    # Initialize an empty list to store all transposed matrices
    site_matrices = []
    for protein_name in config.profile.index.tolist():
        adjusted_pro_site = protein_phos_rate_fitted_analysis(config,protein_name,False)
        site_matrices.append(adjusted_pro_site.drop(config.omics1_name+'_log2FC', axis=1).T)
    phosRate_quantile_matrix = pd.concat(site_matrices,axis=0)
    phosRate_quantile_matrix.to_csv(str(config.outdir)+'/fitted_phosRate_quantile.tsv',sep='\t')
    return phosRate_quantile_matrix

def get_phosRate_quantile_denoised_df(config,phosRate_quantile_matrix):
    frac_para = config.frac_param 
    all_denoised_matrices = []
    for target_site in phosRate_quantile_matrix.index.tolist():
        target_row = phosRate_quantile_matrix.loc[target_site]

        x = np.arange(len(target_row))
        y = target_row.values

        mask = ~ np.isnan(y)
        x_clean = x[mask]
        y_clean = y[mask]
        # # Smooth using Loess; frac controls smoothing level (0-1, larger = smoother)
        loess_smoothed = sm.nonparametric.lowess(y_clean, x_clean, frac=frac_para)  # extract smoothed data
        x_smooth = loess_smoothed[:, 0]
        y_smooth = loess_smoothed[:, 1]
        all_denoised_matrices.append(pd.DataFrame(data=y_smooth).T)
    final_denoised_matrices= pd.concat(all_denoised_matrices,axis=0)
    final_denoised_matrices.index = phosRate_quantile_matrix.index
    final_denoised_matrices.to_csv(str(config.outdir)+'/fitted_phosRate_quantile_denoised.tsv',sep='\t')

def vs_traditional_phos_rate(config,phosRate_quantile_matrix):
    # Calculate classic phosphorylation rate log2FC
    pro_log2FC = config.prodiff['logFC']
    phos_log2FC = config.phosphoprodiff['logFC']

    # Convert phos_log2FC to a DataFrame for merging
    phos_df = phos_log2FC.reset_index()
    phos_df.columns = [config.omics2_name, config.omics2_name+'_log2FC']

    # Convert pro_log2FC to a DataFrame
    pro_df = pro_log2FC.reset_index()
    pro_df.columns = [config.omics1_name, config.omics1_name+'_log2FC']
    # Merge data: first join phosphorylation sites with their corresponding proteins
    merged_df = pd.merge(phos_df, config.phosphopro_pro_cor, on=config.omics2_name, how='left')

    # Then merge in protein log2FC data
    merged_df = pd.merge(merged_df, pro_df, on=config.omics1_name, how='left')
    
    # Compute difference between phosphosite log2FC and protein log2FC
    merged_df['log2FC_'+config.omics2_name] = merged_df[config.omics2_name+'_log2FC'] - merged_df[config.omics1_name+'_log2FC']
    merged_df.dropna(inplace=True)

    site_to_phosrate = merged_df.set_index(config.omics2_name)['log2FC_'+config.omics2_name]
    phosRate_quantile_matrix['total_phos_rate'] = phosRate_quantile_matrix.index.map(site_to_phosrate)
    abs_values = np.abs(phosRate_quantile_matrix.iloc[:, :(config.quantile_rank + 1)])
    phosRate_quantile_matrix['diff_rate'] = (abs_values >= config.phos_rate_FC).sum(axis=1) / (config.quantile_rank + 1)
    mask = (
        (phosRate_quantile_matrix['total_phos_rate'] > -config.phos_rate_FC) &
        (phosRate_quantile_matrix['total_phos_rate'] < config.phos_rate_FC) &
        (phosRate_quantile_matrix['diff_rate'] > config.vs_total_phos_rate_lower) &
        (phosRate_quantile_matrix['diff_rate'] < config.vs_total_phos_rate_upper)
        )
    total_phos_rate_notsig_df = phosRate_quantile_matrix[mask]
    # Sort by diff_rate in descending order
    total_phos_rate_notsig_df = total_phos_rate_notsig_df.sort_values(by='diff_rate', ascending=False)
    total_phos_rate_notsig_df['avergae'] = total_phos_rate_notsig_df.iloc[:, 0:100].mean(axis=1)
    # Merge total_phos_rate_notsig_df with phosphopro_pro_cor to add protein mapping information
    total_phos_rate_notsig_df = pd.merge(total_phos_rate_notsig_df, config.phosphopro_pro_cor, left_index=True, right_on=config.omics2_name, how='left')
    # Restore index to the 'PhosProtein' column (merge removed the original index)
    total_phos_rate_notsig_df.index = total_phos_rate_notsig_df['PhosProtein']
    # Then drop the 'PhosProtein' column from total_phos_rate_notsig_df
    total_phos_rate_notsig_df.drop(columns=['PhosProtein'], inplace=True)

    total_phos_rate_notsig_df.to_csv(str(config.outdir)+'/fitted_phosRate_quantile_vs_total.tsv',sep='\t')

    for top_diff_rate_site in (total_phos_rate_notsig_df.index.tolist()[:config.top_tail_diff_site_num] +  total_phos_rate_notsig_df.index.tolist()[(len(total_phos_rate_notsig_df.index.tolist())-5):len(total_phos_rate_notsig_df.index.tolist())]):
        protein_name_series = config.phosphopro_pro_cor[config.phosphopro_pro_cor[config.omics2_name] == top_diff_rate_site][config.omics1_name]
        total_phos_rate = total_phos_rate_notsig_df.loc[top_diff_rate_site,'total_phos_rate']
        if protein_name_series.empty:
            print(f"Site {top_diff_rate_site} has no corresponding protein, skipping")
            continue
        protein_name = protein_name_series.values[0]  # take the first match
        # print(protein_name+'\t'+top_diff_rate_site)
        config.plot_site_ids = [top_diff_rate_site]
        protein_phos_rate_fitted_analysis(config, protein_name, True,total_phos_rate)

def phos_rate_subtyping(config):
    try:
        # 1. fitting4 and Mapping
        phosRate_quantile_matrix = get_phosRate_quantile_df(config)
        # 2. denoise
        get_phosRate_quantile_denoised_df(config,phosRate_quantile_matrix)
        # 3. clustering [mfuzz]
        def to_r_path(path):
            return str(path).replace("\\", "/")


        # Build cmd — pass all parameters as raw strings, without adding any quotes!
        cmd = [
            "Rscript",
            to_r_path(config.r_script_path),
            "-i", to_r_path(config.outdir / Path('fitted_phosRate_quantile_denoised.tsv')),
            "-o", to_r_path(config.outdir),
            "-c", str(config.cluster_num),
            "-f", str(config.phos_rate_FC),
            "--pluscolor", str(config.Mfuzz_pluscolor),
            "--minuscolor", str(config.Mfuzz_minuscolor),
            "--allcolor", str(config.Mfuzz_allcolor),
            "--onlyupcolor", str(config.Mfuzz_onlyupcolor),
            "--onlydowncolor", str(config.Mfuzz_onlydowncolor),
            "--updowncolor", str(config.Mfuzz_updowncolor),
            "--notsigcolor", str(config.Mfuzz_notsigcolor)
        ]

        full_cmd = cmd

        print(f"📋 Executing command:\n{' '.join(str(arg) for arg in full_cmd)}")
        result = subprocess.run(full_cmd, capture_output=True, text=True)

        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        # 4.vs traditional phos rate
        vs_traditional_phos_rate(config,phosRate_quantile_matrix)
        return True
    except Exception as e:
        print(f"phos_rate_subtyping error: {e}", file=sys.stderr)
        return False



