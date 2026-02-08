#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.process_time_ns()
print('---------------------------------------')
print(time.strftime("RBS_IVcurve.py started: %a, %d %b %Y %H:%M:%S", time.localtime()))
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import h5py
import math
import xraydb
import plotly
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import odrpack as odr
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from getmac import get_mac_address as gma
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
from RBS_functions import *
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
# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)
# plt.rcParams.update({'font.size': 8})

# plt.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     "pgf.preamble": "\n".join([
#           r"\usepackage{mathtools}",
#           r'\usepackage{amsmath}',
#      ]),
# })

plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
color_schemes = load_colors()

end_setup = time.process_time_ns()
elapsed_setup = (end_setup - start_setup)/1e6

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} ms)')
print('---------------------------------------')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def h5_plotter(measurement, m_volt, m_curr, x_err_value, y_err_value, scaling, file_suffix):
    DPI = 300
    plot_type = 'FullData_andFit'
    fig = plt.figure(figsize=(8,4), dpi=DPI, layout='constrained')
    gs = GridSpec(2, 4, figure=fig)
    
    c_index = 0
    x_min, x_max, y_min, y_max = np.min(m_volt), np.max(m_volt), np.min(m_curr), np.max(m_curr)
    # scaling = 10**(-(math.floor(math.log(np.max(m_curr), 10))))
    print(scaling)
    
    #----------------- Measurement Plotter -----------------#
    # for k in list(h5dict.keys()):
    #     measurement = h5dict[k]
    #     scaling = math.floor(math.log(number, 10))
        
    #     if (x_min > np.min(measurement['voltage'])):
    #         x_min = np.min(measurement['voltage'])
            
    #     if (x_max < np.max(measurement['voltage'])):
    #         x_max = np.max(measurement['voltage'])
            
    #     if (y_min > np.min(measurement['current'])):
    #         y_min = np.min(measurement['current'])
        
    #     if (y_max < np.max(measurement['current'])):
    #         y_max = np.max(measurement['current'])
    
    #----------------- Axes Preparer -----------------#
    ax1 = fig.add_subplot(gs[: , :-2]) # I-V plot
    ax2 = fig.add_subplot(gs[: , 2:], sharey=ax1) # I-sqrt(V) plot
    # ax3 = fig.add_subplot(gs[2 , :]) # additional info
    #----------------- Axes Preparer -----------------#
    
    #---------- Measurements ----------# in prefixed ampere (micro/nano)
    ax1.errorbar(x=m_volt,
                y=m_curr*scaling,
                yerr=y_err_value,
                capsize=1,
                capthick=0.4,
                elinewidth=0.2,
                ms=3,
                fmt='x',
                color=color_schemes['c_complementary'][c_index],
                label='measurement')
    
    ax2.errorbar(x=np.sqrt(m_volt),
                y=m_curr*scaling,
                yerr=y_err_value,
                capsize=1,
                capthick=0.4,
                elinewidth=0.2,
                ms=3,
                fmt='x',
                color=color_schemes['c_complementary'][c_index],
                label='measurement')
    
    #---------- Measurements ----------#
    
    #---------- Fitting Procedure ----------#
    p_list =[0.5,10,float(np.min(m_curr))*scaling,0.1]

    b_lower, b_upper = np.array([-10,-10,-10,-10]), np.array([1000,1000,1000,1000])

    
    fit_params = evaluator(sqrt_func=sqrt_func, param_list=p_list, boundary_list=(b_lower,b_upper),
                       x=m_volt, y=m_curr*scaling,
                       xerr=np.array([x_err_value]*len(m_volt)), yerr=np.array([y_err_value]*len(m_curr))) #in prefixed ampere (micro/nano)
    
    fit_params_pureSqrt = [fit_params[0],fit_params[1],0,0]
    print(fit_params_pureSqrt)
    
    # x_lin = np.linspace(np.min(m_volt),np.max(m_volt),5000)
    x_lin = np.linspace(0,np.max(m_volt),5000)
    
    ax1.plot(x_lin,
            sqrt_func(x_lin, fit_params),
            lw=0.75,
            ls='-',
            color='black',
            alpha=0.9,
            zorder=2,
            label='fit'
            )
    
    ax2.plot(np.sqrt(x_lin),
            sqrt_func(x_lin, fit_params),
            lw=0.75,
            ls='-',
            color='black',
            alpha=0.9,
            zorder=2,
            label=r'$f_{\mathrm{fit}}(U_{\mathrm{B}}) = \alpha_1\sqrt{\alpha_2 U_{\mathrm{B}}} + \beta_1 U_{\mathrm{B}} + \beta_2$'
            )
    
    txtout_file_name = f'{measurement['det_type']}_' + plot_type + '_' + f'{measurement['det_id']}' + '_' + f'{file_suffix}.txt'
    
    with open(txtout_file_name,'w') as txtout:
        txtout.write('Fit parameters \n')
        txtout.write(f'alpha_1 = {fit_params[0]}\n')
        txtout.write(f'alpha_2 = {fit_params[1]}\n')
        txtout.write(f'beta_1 = {fit_params[2]}\n')
        txtout.write(f'beta_2 = {fit_params[3]}\n')
    
    #label=r'$f_{\mathrm{fit}}(U_{\mathrm{B}}) = \alpha_1\sqrt{\alpha_2 U_{\mathrm{B}}} + \beta_1 U_{\mathrm{B}} + \beta_2$'
    
    ax2.plot([],[],alpha=0, label=fr'$\alpha_1 = {fit_params[0]*1e3:.3f}$'+ rf'$\cdot 10^{-3}$ / $\alpha_2 = {fit_params[1]:.3f}$')
    ax2.plot([],[],alpha=0, label=fr'$\beta_1 = {fit_params[2]*1e3:.3f}\cdot 10^{-3}$ / $\beta_2 = {fit_params[3]:.3f}$')

    #---------- Fitting Procedure ----------#
    
    c_index += 1
    #----------------- Measurement Plotter -----------------#
    
    #----------------- Detector Image includer -----------------#
    det_pic_file = detector_pic(measurement['det_id'])
    img = plt.imread(det_pic_file)
    annotation = TextArea(f"{measurement['det_type']} Detector \n ID: {measurement['det_id']}", textprops=dict(color="black", fontsize=5, multialignment='center'))
    imagebox = OffsetImage(img, zoom=0.05)
    stacked = VPacker(children=[imagebox, annotation],
                 align="center",
                 pad=0,
                 sep=5)
    
    ab = AnnotationBbox(offsetbox=stacked, xy=(0.15,0.85), xycoords='axes fraction', frameon=True)

    ax1.add_artist(ab)
    #----------------- Detector Image includer -----------------#
    
    
    #----------------- Formulae includer -----------------#
    Formulae_Anno = TextArea(r'$\displaystyle I_{\mathrm{bias}} = \nu\frac{\sqrt{\rho U_{\mathrm{B}}}}{\tau_0}$' + '\n\n' + 
                             r'$\displaystyle I_{\mathrm{bias}}^{\mathrm{Sze}} = \sqrt{\frac{D_p}{\tau_p}} \frac{n_i^2}{N_D} + \nu\frac{\sqrt{\rho U_{\mathrm{B}}}}{\tau_p}$' + '\n\n' + 
                             r'$\displaystyle f_{\mathrm{fit}}(U_{\mathrm{B}}) = \alpha_1\sqrt{\alpha_2 U_{\mathrm{B}}} + \beta_1 U_{\mathrm{B}} + \beta_2$',
                             textprops=dict(fontsize=5, horizontalalignment='center', linespacing=1.5))
    # Fit_Anno = TextArea(fr'$\alpha_1 = {fit_params[0]*1e3:.3f}$' + r'$\cdot 10^{-3}$' + '/' + fr'$\alpha_2 = {fit_params[1]:.3f}$')
    Formulae_stacked = VPacker(children=[Formulae_Anno],
                 align="center",
                 pad=0,
                 sep=5,)
    
    Formulae_ab = AnnotationBbox(offsetbox=Formulae_stacked, xy=(0.19,0.5), xycoords='axes fraction', frameon=True)

    # ax2.add_artist(Formulae_ab)
    #----------------- Formulae includer -----------------#
    
    
    #----------------- Fit Parameter includer -----------------#
    Fit_Anno = TextArea('Fit parameters' + '\n' + fr'$\alpha_1 = {fit_params[0]*1e3:.3f}$' + r'$\cdot 10^{-3}$' + '\n' + fr'$\alpha_2 = {fit_params[1]:.3f}$' + 
                        '\n' +
                        fr'$\beta_1 = {fit_params[2]*1e3:.3f}$' + r'$\cdot 10^{-3}$' + '\n' + fr'$\beta_2 = {fit_params[3]:.3f}$', textprops=dict(fontsize=5))
    # Fit_Anno = TextArea(fr'$\alpha_1 = {fit_params[0]*1e3:.3f}$' + r'$\cdot 10^{-3}$' + '/' + fr'$\alpha_2 = {fit_params[1]:.3f}$')
    Fit_stacked = VPacker(children=[Fit_Anno],
                 align="left",
                 pad=0,
                 sep=5,)
    
    Fit_ab = AnnotationBbox(offsetbox=Fit_stacked, xy=(0.15,0.85), xycoords='axes fraction', frameon=True)

    # ax2.add_artist(Fit_ab)
    #----------------- Fit Parameter includer -----------------#
    
           
    #----------------- Plot Finalization -----------------#
    ax1.annotate('a)',
                 xy=(0.,1.), xycoords='axes fraction',
                 xytext=(+1, -1), textcoords='offset fontsize', fontsize=5)
    ax2.annotate('b)',
                 xy=(0.,1.), xycoords='axes fraction',
                 xytext=(+1, -1), textcoords='offset fontsize', fontsize=5)
    # ax3.annotate('c)',
                #  xy=(0.,1.), xycoords='axes fraction',
                #  xytext=(+1, 0), textcoords='offset fontsize', fontsize=5)
    
    ax1.set_xlabel(r'Bias Voltage / V')
    ax2.set_xlabel(r'Bias Voltage / $\sqrt{\mathrm{V}}$')
    ax1.set_ylabel(r'Leakage Current / $\mu$A')
    # plt.ylabel(r'Leakage Current / nA')
    
    if (x_min < 5 and x_max >= 100):
        x_min = -20
    
    # plt.xlim(x_min*0.5, x_max*1.05)
    ax1.set_xlim(0, x_max*1.05)
    ax2.set_xlim(0, np.sqrt(x_max)*1.05)
    
    ax1.set_ylim((y_min*scaling)*0.99, (y_max*scaling)*1.02)
    
    # plt.ylim((y_min*scaling)*0.5, (y_max*scaling)*1.35)
    # plt.ylim(0.007,4)
    
    # plt.xlim(0,500)
    # plt.ylim(1.8e-6,4e-6)
    # plt.ylim(1.5,4)
    
    plt.yscale('linear')
    
    ax1.grid(which='both')
    ax2.grid(which='both')
    ax1.legend(loc='lower right', fontsize=5)
    ax2.legend(loc='lower right', fontsize=5)
    # plt.tight_layout()
    
    # ax3.xaxis.set_major_formatter(ticker.NullLocator())
    # ax3.xaxis.set_major_formatter(ticker.NullFormatter())
    # ax3.axis("off")
    
    
    save_name_PIIPS = f'{measurement['det_type']}_' + plot_type + '_' + f'{measurement['det_id']}' + '_' + f'{file_suffix}'
    plt.savefig('.//plots//' + save_name_PIIPS + '.pdf', dpi=DPI, transparent=False)
    plt.savefig('.//plots//' + save_name_PIIPS + '.png', dpi=DPI, transparent=False)
    
    plt.show()
    #----------------- Plot Finalization -----------------#
    return 0


if __name__ == "__main__":
    start_routine = time.process_time_ns()
    
    
    # fc = file_collector(dettype='SSB',id='17-440D')
    fc = file_collector(dettype='PIIPS',id='23732')
    h5dict = h5_data_compactor(fc)
    # m_volt, m_curr = h5_measurement_combiner(h5dict)
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        file_suffix = k
        m_volt = measurement['voltage'][1:]
        m_curr = measurement['current'][1:]
        print(m_curr)
        h5_plotter(measurement, m_volt, m_curr, 0.3, 5e-4, 1e6, file_suffix)
    
    print('---------------------------------------')
    end_routine = time.process_time_ns()
    elapsed_routine = (end_routine - start_routine)/1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.4f} s)')
    print('---------------------------------------')