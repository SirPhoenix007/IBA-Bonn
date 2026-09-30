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

from . import polygauss
from . import colors
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

def toy_model(N:int, snr:int, baseline_type:str, plot_flag:bool, D:int=0):
    rng = np.random.default_rng()
    
    x = np.linspace(0,8191,8192)
    # GAUSS
    heights = rng.integers(low=100, high=5e4, size=N)
    widths = rng.integers(low=10, high=100, size=N)
    centers = rng.integers(low=100, high=8100, size=N)
    # print(heights, widths, centers)
    gauss_info = []
    for i in range(0,N):
        gauss_info.append(heights[i])
        gauss_info.append(centers[i])
        gauss_info.append(widths[i])
        
    peak_data = polygauss.multi_gauss(gauss_info, x)
    P_signal = np.mean(peak_data**2)
    
    # BASELINE
    if (baseline_type == 'sine'):
        if (N > 5):
            sine_amp = rng.uniform(low=10000, high=2000*N)
        else:
            sine_amp = rng.uniform(low=10000, high=10000 + 2000*N)
        sine_freq = rng.uniform(low=0.00025, high=0.0006)
        sine_shift = rng.uniform(low=0, high=10000)
        bsl_param = [sine_amp, sine_freq, sine_shift]
        baseline = sine_amp*np.sin(sine_freq*x + sine_shift)
    elif (baseline_type == 'lin'):
        steep = rng.uniform(low=-4, high=4)
        if (steep < 0):
            offset = rng.uniform(low=(np.abs(steep)*8192)+5000, high=(np.abs(steep)*8192)+15000)
        else:
            offset = rng.uniform(low=0, high=10000)
        bsl_param = [offset,steep]
        baseline = offset + x*steep
    elif (baseline_type == 'exp'):
        # a*e^{bx+c}
        exp_a = rng.uniform(low=10, high=20*N)
        exp_b = rng.uniform(low=-1e-3, high=-1e-4)
        exp_c = rng.uniform(low=4, high=8)
        bsl_param = [exp_a, exp_b, exp_c]
        baseline = exp_a*np.exp(exp_b*x + exp_c)
    elif (baseline_type == 'poly'):
        degree = D
        poly_param = rng.uniform(low=-0.2, high=0.2, size=D+1)
        poly = np.polynomial.polynomial.Polynomial(poly_param)
        baseline = poly(x/100)
        bsl_param = poly_param
    
       
    # NOISE
    P_noise = P_signal / (10**(snr/10))
    sigma_noise = np.sqrt(P_noise)
    noise = rng.normal(loc=0.0, scale=sigma_noise, size=8192)
    # print(noise)
    # print(f'SNR: {snr:.3f}')
    
    if (plot_flag == True):
        plt.figure(figsize=(20,4), dpi=250)
        plt.plot(x, peak_data, zorder=4, color='black', label='pure data')
        plt.plot(x, peak_data + noise + baseline, zorder=3, color='firebrick', label='full synthetic data')
        # plt.plot(x, peak_data + baseline + 10000)
        plt.plot(x, baseline, color='teal', zorder=2, label='pure baseline')
        plt.grid(axis='y')
        plt.legend()
        plt.show()
    
    toy_model_data = {
        "Bins": x,
        "NumberOfPeaks": N,
        "SNR": snr,
        "GaussianPeaks": gauss_info,
        "Noise": noise,
        "BaselineType": baseline_type,
        "BaselineParameter": bsl_param,
        "Baseline": baseline,
        "SyntheticData": peak_data+noise+baseline       
    }
    return toy_model_data