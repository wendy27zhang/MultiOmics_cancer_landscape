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
