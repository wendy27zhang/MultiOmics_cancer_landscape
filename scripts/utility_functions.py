### this is functions used for this study

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

### add t_test to the adata object compare left_type with right_type
def t_test(adata_obj, sample_type_col = "SampleType", left_type = "T", right_type="N"):
    left_vector = adata_obj.obs[sample_type_col] == left_type
    right_vector = adata_obj.obs[sample_type_col] == right_type

    mean_left, mean_right, mean_diff, t_vals, p_vals = [], [], [], [], []
    for i in range(adata_obj.X.shape[1]):
        x_left = adata_obj.X[left_vector, i]
        x_right = adata_obj.X[right_vector, i]
        mean_left.append(np.mean(x_left))
        mean_right.append(np.mean(x_right))
        mean_diff.append(np.mean(x_left) - np.mean(x_right))
        t_stat, p_val = ttest_ind(x_left, x_right, equal_var=False)
        t_vals.append(t_stat)
        p_vals.append(p_val)

    new_vars = pd.DataFrame({
        "left_mean" : mean_left,
        "right_mean" : mean_right,
        "diff_mean" : mean_diff,
        "t_vals" : t_vals,
        "p_vals" : p_vals}, index=adata_obj.var.index)
    
    adata_obj.var= adata_obj.var.join(new_vars)  

    return adata_obj


# clean the data frame by require each features has at least col_na_cut non-NA fraction of value, 
#                                   each sample has at least row_na_cut fraction of non-NA value
#   and fill the remain NA with median value across the samples
#   df in samples * features format

def clean_data(X_df, col_na_cut, row_na_cut):
    # Remove features with >30% NA
    col_na_pct = X_df.isna().mean()
    X_df = X_df.loc[:, col_na_pct <= col_na_cut]

    # Remove observations with >30% NA
    row_na_pct = X_df.isna().mean(axis=1)
    X_df = X_df.loc[row_na_pct <= row_na_cut, :]

    # Median imputation
    X_df = X_df.apply(lambda x: x.fillna(x.median()), axis=0)
    return X_df

# example
#df_N = clean_data(df_N, col_na_cut, row_na_cut)
#df_T = clean_data(df_T, col_na_cut, row_na_cut)


### read in the raw data indicatated by cancer type and omic_type
### This function use clean_data to clean data and imput NA 
def read_data(cancer_type, omic_type):
    # Tumor file
    sample_type = "Tumor"
    file_tumor = f"../../MultiOmics_cancer_landscape_data/{cancer_type}/Data/{cancer_type}_{omic_type}_gene_abundance_log2_reference_intensity_normalized_{sample_type}.txt"

    # Normal file
    sample_type = "Normal"
    file_normal = f"../../MultiOmics_cancer_landscape_data/{cancer_type}/Data/{cancer_type}_{omic_type}_gene_abundance_log2_reference_intensity_normalized_{sample_type}.txt"

    # Read data
    df_tumor = pd.read_csv(file_tumor, sep="\t")
    df_normal = pd.read_csv(file_normal, sep="\t")

    # Rename the sample names (orignal sample name are same in tumor and normal file, majority of normal samples are matched normal controls)
    df_normal = df_normal.rename(columns={col: f"{col}_N" for col in df_normal.columns if col != "idx"})
    df_tumor = df_tumor.rename(columns={col: f"{col}_T" for col in df_tumor.columns if col != "idx"})

    # Transposed and clean data, remove 30% NA, and fill the remain with median value in samples
    df_normal.set_index(df_normal.columns[0], inplace=True)
    df_transposed_normal = df_normal.T
    df_transposed_normal_clean = clean_data(X_df= df_transposed_normal, col_na_cut =0.3, row_na_cut =0.3)

    df_tumor.set_index(df_tumor.columns[0], inplace=True)
    df_transposed_tumor = df_tumor.T
    df_transposed_tumor_clean = clean_data(X_df= df_transposed_tumor, col_na_cut =0.3, row_na_cut =0.3)

    df_transposed_merged = pd.concat([df_transposed_normal_clean, df_transposed_tumor_clean], axis=0, join="inner")

    return df_transposed_merged
