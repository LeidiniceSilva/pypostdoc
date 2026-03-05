# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "March 02, 2026"
__description__ = "This script compute trend"

import numpy as np
from scipy.stats import t as student_t


def calculate_trend(data, alpha=0.05, per_decade=False):

    nt, ny, nx = data.shape
    t = np.arange(nt, dtype=float)
    
    trend = np.full((ny, nx), np.nan)
    sig_mask = np.full((ny, nx), False)
    
    for i in range(ny):
        for j in range(nx):
            y = data[:, i, j]
            valid = np.isfinite(y)
            
            if np.sum(valid) >= 10:  
                t_valid = t[valid]
                y_valid = y[valid]
                
                # Linear regression
                A = np.vstack([t_valid, np.ones_like(t_valid)]).T
                slope, intercept = np.linalg.lstsq(A, y_valid, rcond=None)[0]
                
                # Calculate trend 
                trend_val = slope * 10.0 if per_decade else slope
                trend[i, j] = trend_val
                
                # Calculate p-value
                y_pred = slope * t_valid + intercept
                residuals = y_valid - y_pred
                ss_res = np.sum(residuals**2)
                
                if len(t_valid) > 2 and ss_res > 0:
                    ss_tot = np.sum((y_valid - np.mean(y_valid))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    df = len(t_valid) - 2
                    
                    # t-statistic
                    t_stat = slope / np.sqrt(ss_res / (df * np.sum((t_valid - np.mean(t_valid))**2)))
                    p_value = 2 * (1 - student_t.cdf(np.abs(t_stat), df))
                    
                    if p_value < alpha:
                        sig_mask[i, j] = True
    
    return trend, sig_mask


