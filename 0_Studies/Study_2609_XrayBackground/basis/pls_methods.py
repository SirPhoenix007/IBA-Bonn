#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import tqdm
import warnings
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import multiprocessing as mup
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from pybaselines import Baseline
from sklearn.metrics import root_mean_squared_error

from . import *
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times"],
    "text.usetex": True,
    "font.size": 8,
    "pgf.rcfonts": False
})


plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
color_schemes = colors.load_colors()

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def arpls_baseline(bin_data:list, bins:list, lam:int = 1e2):
    bsl_fitter = Baseline(x_data=bins)

    try:
        baseline, params = bsl_fitter.arpls(bin_data, lam=lam, max_iter=500)
    except np.linalg.LinAlgError:
        # Extremely large lam can make the banded Whittaker system
        # numerically non-positive-definite (floating point round-off).
        # Treat this lam as a failed fit rather than crashing the worker.
        return None, None
    subtracted = bin_data - baseline
    return baseline, subtracted

def aspls_baseline(bin_data:list, bins:list, lam:int = 1e2):
    bsl_fitter = Baseline(x_data=bins)

    try:
        baseline, params = bsl_fitter.aspls(bin_data, lam=lam, max_iter=500)
    except np.linalg.LinAlgError:
        return None, None
    subtracted = bin_data - baseline
    return baseline, subtracted

def evaluate_baseline(toy_model_data:dict):
    
    plt.figure(figsize=(20,4), dpi=250)
    plt.grid(axis='y', which='both')
    
    lambda_range = np.linspace(2,13,2201)
    
    bins = toy_model_data['Bins']
    data = toy_model_data['SyntheticData']
    real_bsl = toy_model_data['Baseline']
    lam_min = [0,0]
    minimum = [1e5,1e5]
    res_list = [[],[]]
    for lam in tqdm.tqdm(lambda_range):
        arpls_bsl,_ = arpls_baseline(data, bins, 10**lam)
        if arpls_bsl is None:
            res_list[0].append(np.nan)
        else:
            arpls_result = fit_quality.rmse(real_bsl,arpls_bsl, 8192)
            res_list[0].append(arpls_result)
            if arpls_result < minimum[0]:
                minimum[0] = arpls_result
                lam_min[0] = lam
        
        aspls_bsl,_ = aspls_baseline(data, bins, 10**lam)
        if aspls_bsl is None:
            res_list[1].append(np.nan)
        else:
            aspls_result = fit_quality.rmse(real_bsl,aspls_bsl, 8192)
            res_list[1].append(aspls_result)
            if aspls_result < minimum[1]:
                minimum[1] = aspls_result
                lam_min[1] = lam
            
        # print(f'{lam:.4f}: {result:.4f}')
        # plt.plot(bins,arpls_bsl)
    # plt.plot(bins,real_bsl, color='black')
    plt.plot(lambda_range,res_list[0],label='arPLS')
    plt.plot(lambda_range,res_list[1],label='asPLS')
    plt.legend()
    plt.ylim(1,20000)
    plt.yscale('log')
    plt.show()
    # print(f'arPLS: Minimum: {lam_min[0]:.2f}: {minimum[0]:.2f} / {(minimum[0]/(data.mean())):.3f}')
    # print(f'asPLS: Minimum: {lam_min[1]:.2f}: {minimum[1]:.2f} / {(minimum[1]/(data.mean())):.3f}')
    
    
    # as_min,_ = aspls_baseline(data, bins, 10**lam_min[1])
    # ar_min,_ = arpls_baseline(data, bins, 10**lam_min[0])
    
    # plt.figure(figsize=(20,4), dpi=250)
    # plt.grid(axis='y', which='both')
    # plt.plot(bins, real_bsl, label='real Baseline')
    # plt.plot(bins, as_min, label='asPLS Baseline')
    # plt.plot(bins, ar_min, label='arPLS Baseline')
    # plt.legend()
    return {
        "arpls_lambda_min": round(lam_min[0],3),
        "arpls_rmse_min": minimum[0],
        "arpls_rmse_all": res_list[0],
        "aspls_lambda_min": round(lam_min[1],3),
        "aspls_rmse_min": minimum[1],
        "aspls_rmse_all": res_list[1],}