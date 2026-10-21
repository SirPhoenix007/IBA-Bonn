lazy import numpy as np
lazy import scipy as sc
lazy import pandas as pd
lazy import json
lazy import matplotlib.pyplot as plt
lazy import itertools as ito

from iba_bonn_pixe import *

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],3)})
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


def read_ecal():
    with open('./energy_peaks/peaks_for_calibration.json', 'r') as f:
        data = json.load(f)
    print(f'INFO: ENERGY CALIBRATION DATA LOADED FROM FILE')
    return data

def peak_fitter(order:int,lines:list, peaks:list, peaks_err:list, init_values:list):
    
    lines_err = np.array([1]*len(lines))
    
    if order == 1:
        beta = pixe_functions.evaluator_scipy(func=pixe_functions.lin_func, beta0_list=init_values,
                                   x=peaks, y=lines, xerr=peaks_err, yerr=lines_err)
    elif order == 2:    
        beta = pixe_functions.evaluator_scipy(func=pixe_functions.energy_func, beta0_list=init_values,
                           x=peaks, y=lines, xerr=peaks_err, yerr=lines_err)
    
    return beta

def lowerLim_lin(beta, x):
    param = beta[0]
    errors = beta[1]
    return (param[0]-errors[0]) + (param[1]-errors[1])*x

def lowerLim_poly2(beta, x):
    param = beta[0]
    errors = beta[1]
    return (param[0]-errors[0]) + (param[1]-errors[1])*x + (param[2]-errors[2])*x**2

def upperLim_lin(beta, x):
    param = beta[0]
    errors = beta[1]
    return (param[0]+errors[0]) + (param[1]+errors[1])*x

def upperLim_poly2(beta, x):
    param = beta[0]
    errors = beta[1]
    return (param[0]+errors[0]) + (param[1]+errors[1])*x + (param[2]+errors[2])*x**2

def parameter_readback(params:list, order:int):
    if order == 1:
        print('--------------------------------------')
        print('LINEAR FIT:')
        print(f'({round(params[0][0],3)} \u00B1 {round(params[1][0],3)}) eV / {np.abs(round(100*params[1][0]/params[0][0],2))} %')
        print(f'({round(params[0][1],3)} \u00B1 {round(params[1][1],3)}) eV/Nc / {np.abs(round(100*params[1][1]/params[0][1],2))} %')
        print('--------------------------------------')
    elif order == 2:
        print('--------------------------------------')
        print('POLY-2 FIT:')
        print(f'({round(params[0][0],3)} \u00B1 {round(params[1][0],3)}) eV / {np.abs(round(100*params[1][0]/params[0][0],2))} %')
        print(f'({round(params[0][1],3)} \u00B1 {round(params[1][1],3)}) eV/Nc / {np.abs(round(100*params[1][1]/params[0][1],2))} %')
        print(f'({round(params[0][2]*1e6,3)} \u00B1 {round(params[1][2]*1e6,3)}) 10^-6 eV²/Nc² / {np.abs(round(100*params[1][2]/params[0][2],2))} %')
        print('--------------------------------------')

def chi_readback(lines:list, mca_v:list, params:list, order:int):
    bin_data = lines
    if order == 1:
        fit_data = pixe_functions.lin_func(params, np.array(mca_v))
        chis = fit_quality.both_chi_squared(bin_data, fit_data, params)
        print('--------------------------------------')
        print(f'\u03C7² = {chis[0]}')
        print(f'reduced \u03C7² = {chis[1]}')
    elif order == 2:
        fit_data = pixe_functions.energy_func(params, np.array(mca_v))
        chis = fit_quality.both_chi_squared(bin_data, fit_data, params)
        print('--------------------------------------')
        print(f'\u03C7² = {chis[0]}')
        print(f'reduced \u03C7² = {chis[1]}')

def single_setting_ecal(MCA:int, EG:int):
    
    DPI = 250
    
    pfc = read_ecal()
    ecal_data = pfc[str(MCA)][str(EG)]
    ecal_keys = list(ecal_data.keys())
    n = int(len(ecal_data)/3)
    
    try:
        len(ecal_data)%3 == 0
    except:
        raise ValueError(f'The Ecal-data is missing at least one list contribution for {MCA} / {EG}.')
    
    lines, mca_values, mca_errors = [],[],[]
    
    for i in range(0,n):
        lines.append(ecal_data[ecal_keys[i]])
        mca_values.append(ecal_data[ecal_keys[i+n]])
        mca_errors.append(ecal_data[ecal_keys[i+(2*n)]])
    
    lines_flat = list(ito.chain.from_iterable(lines))
    mca_v_flat = list(ito.chain.from_iterable(mca_values))
    mca_e_flat = list(ito.chain.from_iterable(mca_errors))
    print(lines_flat)
    
    init_values = [100,8,0]
    beta_lin= peak_fitter(order=1, lines=lines_flat, peaks=mca_v_flat, peaks_err=mca_e_flat, init_values=init_values)
    beta_poly2= peak_fitter(order=2, lines=lines_flat, peaks=mca_v_flat, peaks_err=mca_e_flat, init_values=init_values)

    beta_export_lin = [beta_lin['param'],beta_lin['errors']]
    beta_export_poly2 = [beta_poly2['param'],beta_poly2['errors']]
    
    parameter_readback(beta_export_lin, order=1)
    chi_readback(lines_flat, mca_v_flat, beta_export_lin[0], order=1)
    parameter_readback(beta_export_poly2, order=2)
    chi_readback(lines_flat, mca_v_flat, beta_export_poly2[0], order=2)
    
    print('ENERGY CALIBRATION: Fit done.')
    
    fig, ax = plt.subplots(figsize=(5,3), dpi=DPI)
    ax.scatter(mca_v_flat, lines_flat, s=10, marker='s', lw=0.1, color=color_schemes['c_rainbow'][6], zorder=3)
    
    xlin = np.linspace(0,8192,8192)
    ax.plot(xlin, pixe_functions.lin_func(beta_export_lin[0],xlin), color='red', lw=0.75, label=r'linear fit', zorder=2)
    ax.fill_between(xlin,lowerLim_lin(beta_export_lin,xlin),upperLim_lin(beta_export_lin,xlin),alpha=0.25, color='red', zorder=2, label=r'$1\sigma$-error')
    ax.plot(xlin, pixe_functions.energy_func(beta_export_poly2[0],xlin), color='green', lw=0.75, label=r'$2^{\mathrm{nd}}$ order polynomial fit', zorder=2)
    ax.fill_between(xlin,lowerLim_poly2(beta_export_poly2,xlin),upperLim_poly2(beta_export_poly2,xlin),alpha=0.25, color='green', zorder=2, label=r'$1\sigma$-error')
        
    ax.set_xlabel(r'MCA channel $N_{\mathrm{c}}$')
    ax.set_ylabel('Characteristic X-ray energy in keV')
    
    ax.set_xlim(0,8193)
    ax.set_xticks(np.arange(0,8193,1024),np.arange(0,8193,1024))
    ax.set_ylim(0,20000)
    ax.set_yticks([0,5000,10000,15000,20000],[0,5,10,15,20])
    
    plt.show()
    
    return 10


single_setting_ecal(8192,3400)

