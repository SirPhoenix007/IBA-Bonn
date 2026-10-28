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

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from basis import *
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
pls_toymodel_location = "C:/Users/schum/Documents/Filing Cabinet/5_PLSSimulationFiles/simulation"
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def read_json_formatted_file(filepath, encoding="utf-8"):
    try:
        return json.load(open(filepath, "r", encoding=encoding))
    except json.JSONDecodeError as e:
        raise ValueError(f"File content is not valid JSON: {e}") from e

'''
simulation_data[i][j][ keys ]

- i: number of experiments: = len(simulation_data) = 100 for most
- j: 
  - 0: toy model data
  - 1: PLS method data
- keys:
  - if j = 0: 'Bins', 'NumberOfPeaks', 'SNR', 'GaussianPeaks', 'Noise', 'BaselineType', 'BaselineParameter', 'Baseline', 'SyntheticData'
  - if j = 1: 'arpls_lambda', 'arpls_rmse', 'aspls_lambda', 'aspls_rmse' (each optimal parameters for given toy model)
'''


# data_json = pls_toymodel_location + "./toy_model_data/toy_model_run_170626_190259.json"

# test_data = read_json_formatted_file(data_json)[0][0]


def plot_uncompressed_PLS(test_data, eval_result):
    print('arPLS:' + f'{eval_result['arpls_lambda_min']}: {eval_result['arpls_rmse_min']}')
    print('asPLS:' + f'{eval_result['aspls_lambda_min']}: {eval_result['aspls_rmse_min']}')

    arPLS,_ = pls_methods.arpls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=10**eval_result['arpls_lambda_min'])
    asPLS,_ = pls_methods.aspls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=10**eval_result['aspls_lambda_min'])
    plt.figure(figsize=(4,2), dpi=300)
    plt.plot(test_data['Bins'],test_data['SyntheticData'], lw=0.5, label='synthetic spectrum')
    plt.plot(test_data['Bins'],arPLS, lw=0.5, label='arPLS')
    plt.plot(test_data['Bins'],asPLS, lw=0.5, label='asPLS')
    # plt.plot(test_data['Bins'],delta)

    # plt.ylim(1.5,2.8)
    plt.grid()
    plt.yscale('linear')
    plt.legend()
    plt.show()
    return 10


def plot_compressed_PLS(comp_test_data, comp_eval_result):
    print('arPLS:' + f'{comp_eval_result['arpls_lambda_min']}: {comp_eval_result['arpls_rmse_min']}')
    print('asPLS:' + f'{comp_eval_result['aspls_lambda_min']}: {comp_eval_result['aspls_rmse_min']}')

    arPLS,_ = pls_methods.arpls_baseline(comp_test_data['SyntheticData'],comp_test_data['Bins'],lam=10**comp_eval_result['arpls_lambda_min'])
    asPLS,_ = pls_methods.aspls_baseline(comp_test_data['SyntheticData'],comp_test_data['Bins'],lam=10**comp_eval_result['aspls_lambda_min'])

    plt.figure(figsize=(4,2), dpi=150)
    plt.plot(comp_test_data['Bins'],comp_test_data['SyntheticData'], lw=0.5, label='compressed synth. spectrum')
    plt.plot(comp_test_data['Bins'],arPLS, lw=0.5, label='arPLS')
    plt.plot(comp_test_data['Bins'],asPLS, lw=0.5, label='asPLS')
    # plt.plot(test_data['Bins'],delta)
    plt.legend()
    # plt.ylim(1.5,2.8)
    plt.yscale('linear')
    plt.grid()
    plt.show()
    return 20


def LL_compression_report(test_data):
    original = test_data['SyntheticData']
    
    comp = compression.LL_Compression(test_data['SyntheticData'])
    base = compression.LL_Compression(test_data['Baseline'])
    
    test_data['SyntheticData'] = comp
    test_data['Baseline'] = base
    decomp = compression.LL_Decompression(comp)
    
    delta_list = original-decomp
    delta = np.sum((original - decomp)/original)
    delta_rmse = root_mean_squared_error(original, decomp)
    
    print('-------------------------------')
    print(f'LL compression/decompression delta: {delta*100:.3g} % / RMSE: {delta_rmse:.3g}')
    print('-------------------------------')
    
    comp_eval_result = pls_methods.evaluate_baseline(test_data, compressed_flag=True)
    
    plot_compressed_PLS(test_data, comp_eval_result)
    return 30
    
def LLS_compression_report(test_data):
    original = test_data['SyntheticData']
    
    comp = compression.LLS_Compression(test_data['SyntheticData'])
    base = compression.LLS_Compression(test_data['Baseline'])
    
    test_data['SyntheticData'] = comp
    test_data['Baseline'] = base
    decomp = compression.LLS_Decompression(comp)
    
    delta_list = original-decomp
    delta = np.sum((original - decomp)/original)
    delta_rmse = root_mean_squared_error(original, decomp)
    
    print('-------------------------------')
    print(f'LLS compression/decompression delta: {delta*100:.3g} % / RMSE: {delta_rmse:.3g}')
    print('-------------------------------')
    
    comp_eval_result = pls_methods.evaluate_baseline(test_data, compressed_flag=True)
    
    plot_compressed_PLS(test_data, comp_eval_result)
    return 40


if __name__ == "__main__":
    test_data = generate_toymodel.toy_model(N=12, snr=25, baseline_type='lin', plot_flag=True)
    eval_result = pls_methods.evaluate_baseline(test_data)
    
    plot_uncompressed_PLS(test_data, eval_result)
    
    LL_compression_report(test_data)
    LLS_compression_report(test_data)
    