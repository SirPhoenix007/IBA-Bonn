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

test_data = generate_toymodel.toy_model(N=12, snr=25, baseline_type='lin', plot_flag=True)



eval_result = pls_methods.evaluate_baseline(test_data)

print('arPLS:' + f'{eval_result['arpls_lambda_min']}: {eval_result['arpls_rmse_min']}')
print('asPLS:' + f'{eval_result['aspls_lambda_min']}: {eval_result['aspls_rmse_min']}')

arPLS,_ = pls_methods.arpls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=eval_result['arpls_lambda_min'])
asPLS,_ = pls_methods.aspls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=eval_result['aspls_lambda_min'])
plt.figure(figsize=(8,4), dpi=300)
plt.plot(test_data['Bins'],test_data['SyntheticData'])
plt.plot(test_data['Bins'],arPLS)
plt.plot(test_data['Bins'],asPLS)
# plt.plot(test_data['Bins'],delta)

# plt.ylim(1.5,2.8)
plt.yscale('linear')

plt.show()

comp = ddlc.Dynamic_Double_Logarithmic_Compression(test_data['SyntheticData'])
base = ddlc.Dynamic_Double_Logarithmic_Compression(test_data['Baseline'])
test_data['SyntheticData'] = comp
test_data['Baseline'] = base
decomp = ddlc.Dynamic_Double_Logarithmic_Decompression(comp)
delta = test_data['SyntheticData'] - decomp



eval_result = pls_methods.evaluate_baseline(test_data)

print('arPLS:' + f'{eval_result['arpls_lambda_min']}: {eval_result['arpls_rmse_min']}')
print('asPLS:' + f'{eval_result['aspls_lambda_min']}: {eval_result['aspls_rmse_min']}')

arPLS,_ = pls_methods.arpls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=eval_result['arpls_lambda_min'])
asPLS,_ = pls_methods.aspls_baseline(test_data['SyntheticData'],test_data['Bins'],lam=eval_result['aspls_lambda_min'])

plt.figure(figsize=(8,4), dpi=300)
plt.plot(test_data['Bins'],comp)
plt.plot(test_data['Bins'],arPLS)
plt.plot(test_data['Bins'],asPLS)
# plt.plot(test_data['Bins'],delta)

# plt.ylim(1.5,2.8)
plt.yscale('linear')

plt.show()