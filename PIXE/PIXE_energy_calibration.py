lazy import numpy as np
lazy import scipy as sc
lazy import pandas as pd
lazy import json
lazy import math
lazy import matplotlib.pyplot as plt
lazy import matplotlib.gridspec as gridspec
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
        return chis
    elif order == 2:
        fit_data = pixe_functions.energy_func(params, np.array(mca_v))
        chis = fit_quality.both_chi_squared(bin_data, fit_data, params)
        print('--------------------------------------')
        print(f'\u03C7² = {chis[0]}')
        print(f'reduced \u03C7² = {chis[1]}')
        return chis

def single_setting_ecal(MCA:int, EG:int):
    
    DPI = 250
    
    xlin = np.linspace(0,MCA,MCA)
    
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
    chis_lin = chi_readback(lines_flat, mca_v_flat, beta_export_lin[0], order=1)
    parameter_readback(beta_export_poly2, order=2)
    chis_poly2 = chi_readback(lines_flat, mca_v_flat, beta_export_poly2[0], order=2)
    
    print('ENERGY CALIBRATION: Fit done.')
    
    fig = plt.figure(figsize=(7, 3), constrained_layout=True, dpi=DPI)
    gs = fig.add_gridspec(3, 3)
    ax = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2, :])
    #TODO: Add Errors to values in plot.
    ax.scatter(mca_v_flat, lines_flat, s=10, marker='s', lw=0.1, color='black', zorder=3)
    
    ax.plot(xlin, pixe_functions.lin_func(beta_export_lin[0],xlin), color=color_schemes['c_rainbow'][1], lw=0.75, label='linear fit:\n ' + f'offset = {round(beta_export_lin[0][0],2)} eV / gradient = {round(beta_export_lin[0][1],3)}' + r' eV/$N_{\mathrm{c}}$', zorder=2)
    ax.fill_between(xlin,lowerLim_lin(beta_export_lin,xlin),upperLim_lin(beta_export_lin,xlin),alpha=0.25, color=color_schemes['c_rainbow'][1], zorder=2, label=r'$1\sigma$-error')
    
    ax.plot(xlin, pixe_functions.energy_func(beta_export_poly2[0],xlin), color=color_schemes['c_rainbow'][6], lw=0.75, label='quadratic fit:\n ' + f'offset = {round(beta_export_poly2[0][0],2)} eV / gradient = {round(beta_export_poly2[0][1],3)}' + r' eV/$N_{\mathrm{c}}$ ' + f' / non-linear =  {round(beta_export_poly2[0][2]*10**6,3)}' + r'$\cdot 10^{-6}$ eV$^2$/$N_{\mathrm{c}}^2$', zorder=2)
    ax.fill_between(xlin,lowerLim_poly2(beta_export_poly2,xlin),upperLim_poly2(beta_export_poly2,xlin),alpha=0.25, color=color_schemes['c_rainbow'][6], zorder=2, label=r'$1\sigma$-error')
    
    # ax.plot([],[],label=f'{round(beta_export_lin[0][0],2)} + {round(beta_export_lin[0][1],3)}' + r'$N_{\mathrm{c}}$/eV')
    
    #--------------------------------------------------------------
    
    res_lin = pixe_functions.lin_func(beta_export_lin[0],np.array(mca_v_flat))-np.array(lines_flat)
    res_poly2 = pixe_functions.lin_func(beta_export_poly2[0],np.array(mca_v_flat))-np.array(lines_flat)
    
    ax2.scatter(mca_v_flat, res_lin, s=10, marker='s', lw=0.1, color=color_schemes['c_rainbow'][1], zorder=3, label=r'residuals from linear fit: $\chi^2$ = ' + f'{chis_lin[0]:.4g}')
    ax2.vlines(mca_v_flat, 0, res_lin, color=color_schemes['c_rainbow'][1], alpha=0.7, lw=0.4)
    
    ax2.scatter(mca_v_flat, res_poly2, s=10, marker='o', lw=0.1, color=color_schemes['c_rainbow'][6], zorder=3, label=r'residuals from quadratic fit: $\chi^2$ = ' + f'{chis_poly2[0]:.4g}')
    ax2.vlines(mca_v_flat, 0, res_poly2, color=color_schemes['c_rainbow'][6], alpha=0.7, lw=0.4)
    ax2.axhline(0, color="black", linewidth=1)
        
    # ax.set_xlabel(r'MCA channel $N_{\mathrm{c}}$')
    ax.set_ylabel('Characteristic X-ray energy in keV')
    
    ax2.set_xlabel(r'MCA channel $N_{\mathrm{c}}$')
    ax2.set_ylabel('Residuals')
    
    ax.set_xlim(0,MCA+1)
    ax.set_xticks(np.arange(0,MCA+1,1024),np.arange(0,MCA+1,1024))
    
    ax2.set_xlim(0,MCA+1)
    ax2.set_xticks(np.arange(0,MCA+1,1024),np.arange(0,MCA+1,1024))
    
    
    ymax = math.ceil(np.max(lines_flat) / 5000) * 5000 + 5000 
    ax.set_ylim(0,ymax - 2500)
    ax.set_yticks(np.arange(0,ymax,5000),np.arange(0,ymax/1000,5))
    
    resi_ymax = math.ceil(max(np.max(res_lin),np.max(res_poly2))/50)*50
    resi_ymin = math.floor(min(np.min(res_lin),np.min(res_poly2))/50)*50
    # ax2.set_ylim(resi_ymin*1.1,resi_ymax*1.1)
    ax2.set_yticks(np.arange(resi_ymin,resi_ymax, 50))
    
    ax.grid(visible=True)
    ax2.grid(visible=True)
    
    ax.legend(fontsize=6)
    ax2.legend(fontsize=6)
    
    plt.show()
    
    return 10


single_setting_ecal(8192,5461) #TODO: Create second function with identical plot but for slides

