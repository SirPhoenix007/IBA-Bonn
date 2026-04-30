
#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.process_time_ns()
print('---------------------------------------')
print(time.strftime("PIXE_plotting.py started: %a, %d %b %Y %H:%M:%S", time.localtime()))
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import h5py
import math
import tqdm
import xraydb
import plotly
import argparse
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
from scipy.special import voigt_profile
from getmac import get_mac_address as gma
from itertools import chain
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
from PIXE_functions import *
from PIXE_polygauss import *
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
color_schemes = load_colors()


end_setup = time.process_time_ns()
elapsed_setup = (end_setup - start_setup)/1e6

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} ms)')
print('---------------------------------------')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
# PARSER SETUP
parser = argparse.ArgumentParser(description='PIXE Spectrum Plotting')

parser.add_argument('-sm','--simple_mode', choices=['None','single', 'dual'], default='None', help='Plotting mode for simple plots: single or dual axis')
parser.add_argument('-dm','--detailed_mode', choices=['None','1','2'], default='None', help='Plotting mode for detailed plots: single axis + expected lines (1,2)')

parser.add_argument('-l','--log', help='Plot y-axis in logarithmic scale', action='store_true')

parser.add_argument('-p','--path', help='Directory of data to insert into routine.')
parser.add_argument('-s','--save', help='Save the generated plot as .png and .pdf files', action='store_true')


#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~#

def read_ecal(ecal_data:str):
    with open(ecal_data, 'r') as f:
        data = json.load(f)
    print(f'INFO: ENERGY CALIBRATION DATA LOADED FROM {ecal_data}')
    return data

def read_full_directory(directory_name:str):
    f_c = file_collector(directory_name)
    return f_c

def measurement_parameters(file:str, ecal_data:dict):
    curr_data_set = read_json_formatted_file(file, encoding="utf-8")
    curr_raw_data = curr_data_set['RawData']
    curr_config = curr_data_set['Configuration']
    curr_parameters = curr_config['Parameters']
    curr_MCA = int(curr_parameters[2]['ParameterValue'])
    curr_EG = int(curr_parameters[16]['ParameterValue'])
    curr_file_name = file.split('//')[-1]
    print(f'INFO: FILE {file} READ (MCA: {curr_MCA}, EG: {curr_EG})')
    
    for i in range(1,6):
        if (curr_MCA == ecal_data[str(i)]['MCA'] or 2*curr_MCA == ecal_data[str(i)]['MCA']) and curr_EG == ecal_data[str(i)]['EG']:
            ecal_MCA = ecal_data[str(i)]['MCA']
            ecal_EG = ecal_data[str(i)]['EG']
            ecal_param = ecal_data[str(i)]['param']
            ecal_paramUn = ecal_data[str(i)]['paramUn']
            print(f'INFO: ECAL for FILE: DC-offset {ecal_param[0]:.4f} | Linear {ecal_param[1]:.4f} | Quadratic {ecal_param[2]*1e6:.4f} x 10^-6')
            continue
    return curr_file_name, curr_raw_data, curr_MCA, curr_EG, ecal_param, ecal_paramUn

def data_converter(measurement_data:list, ecal_param:list):
    converted_data = []
    for i in range(len(measurement_data)):
        converted_value = ecal_param[0] + ecal_param[1]*measurement_data[i] + ecal_param[2]*measurement_data[i]**2
        converted_data.append(converted_value)
    return converted_data

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~#

def plot_basic_singleAxis_spectrum(meas_file:str, meas_EG:int, meas_MCA:int, meas_data:list, ecal_param:list, save_flag:bool=False, log_flag:bool=False):
    '''
    Single axis plot with energy vs. counts.\n
    Simplified xticks for readable energy values.\n 
    No secondary x-axis with MCA channels.
    '''
       
    DPI = 200
    
    if (len(meas_data) == 8192):
        xticks_list = np.arange(0,8193,1024)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 4096):
        xticks_list = np.arange(0,4097,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 2048):
        xticks_list = np.arange(0,2049,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
        
    if convxticks_list[-1] <= 10000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/500) *500, 500)
    elif convxticks_list[-1] <= 25000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/1000) *1000, 1000)
    elif convxticks_list[-1] <= 50000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/2500) *2500, 2500)
    elif convxticks_list[-1] <= 100000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/5000) *5000, 5000)    
        
    meas_file_name = meas_file.split('.')[0]
    meas_data = meas_data[:-1]
    
    converted_x_axis = data_converter(np.linspace(0, meas_MCA-2, meas_MCA-1), ecal_param)
    
    fig, ax = plt.subplots(figsize=(6,4), dpi=DPI)
    ax.plot(converted_x_axis, meas_data, color=color_schemes['c_five2'][0], label=f'Calibrated Spectrum \n ID: {meas_file_name}')
    
    ax.set_xlim(0, convxticks_list[-1])
    ax.set_xticks(eV_xticks_list, [f'{int(i/1000)}' for i in eV_xticks_list])
    
    if log_flag == True:
        ax.set_yscale('log')
        ax.set_ylim(1e0,np.max(meas_data)*1.2)
    else:
        ax.set_yscale('linear')
        ax.set_ylim(0,np.max(meas_data)*1.2)
    
    plt.xlabel('Energy / keV', fontsize=10)
    plt.ylabel('Counts', fontsize=10)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    
    if (save_flag == True and log_flag == True):
        plt.savefig(f'./plots/calibrated/only_data_log/{meas_file_name}_singleAx_log.png', transparent=False, dpi=DPI)
        plt.savefig(f'./plots/calibrated/only_data_log/{meas_file_name}_singleAx_log.pdf', transparent=False, dpi=DPI)
    elif (save_flag == True and log_flag == False):
        plt.savefig(f'./plots/calibrated/only_data_lin/{meas_file_name}_singleAx.png', transparent=False, dpi=DPI)
        plt.savefig(f'./plots/calibrated/only_data_lin/{meas_file_name}_singleAx.pdf', transparent=False, dpi=DPI)
    else:
        plt.show()

    plt.close()
    return 0

def plot_basic_dualAxis_spectrum(meas_file:str, meas_EG:int, meas_MCA:int, meas_data:list, ecal_param:list, save_flag:bool=False, log_flag:bool=False):
    '''
    Dual axis plot with energy vs. counts.\n
    xticks depending on MCA channels.\n 
    Secondary x-axis with MCA channels.
    '''
    
    DPI = 200
    
    if (len(meas_data) == 8192):
        xticks_list = np.arange(0,8193,1024)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 4096):
        xticks_list = np.arange(0,4097,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 2048):
        xticks_list = np.arange(0,2049,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
    
    
    meas_file_name = meas_file.split('.')[0]
    meas_data = meas_data[:-1]
    
    converted_x_axis = data_converter(np.linspace(0, meas_MCA-2, meas_MCA-1), ecal_param)
    
    fig, ax = plt.subplots(figsize=(6,4), dpi=DPI)
    ax.plot(converted_x_axis, meas_data, color=color_schemes['c_five2'][0], label=f'Calibrated Spectrum \n ID: {meas_file_name}')
    
    ax.set_xlim(0, convxticks_list[-1])
    ax.set_xticks(convxticks_list, [f'{int(i)}' for i in convxticks_list])
    
    if log_flag == True:
        ax.set_yscale('log')
        ax.set_ylim(1e0,np.max(meas_data)*1.2)
    else:
        ax.set_yscale('linear')
        ax.set_ylim(0,np.max(meas_data)*1.2)
    
    secax = ax.secondary_xaxis('top')
    secax.set_xlim(0,xticks_list[-1])
    secax.set_xticks(convxticks_list,xticks_list)
    secax.set_xlabel('MCA Channel', fontsize=10)
    
    plt.xlabel('Energy / eV', fontsize=10)
    plt.ylabel('Counts', fontsize=10)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    
    if (save_flag == True and log_flag == True):
        plt.savefig(f'./plots/calibrated/only_data_log/{meas_file_name}_singleAx_log.png', transparent=False, dpi=DPI)
        plt.savefig(f'./plots/calibrated/only_data_log/{meas_file_name}_singleAx_log.pdf', transparent=False, dpi=DPI)
    elif (save_flag == True and log_flag == False):
        plt.savefig(f'./plots/calibrated/only_data_lin/{meas_file_name}_singleAx.png', transparent=False, dpi=DPI)
        plt.savefig(f'./plots/calibrated/only_data_lin/{meas_file_name}_singleAx.pdf', transparent=False, dpi=DPI)
    else:
        plt.show()
    
    plt.close()
    return 0

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~#

def plot_data_expected_singleAxis_spectrum(meas_file:str, meas_MCA:int, meas_data:list, ecal_param:list, elem_symbol:str, both_axes:bool=False, save_flag:bool=False, log_flag:bool=False):
    '''
    Single axis plot with energy vs. counts.\n
    Simplified xticks for readable energy values.\n 
    No secondary x-axis with MCA channels.
    '''
       
    DPI = 200
    
    if (len(meas_data) == 8192):
        xticks_list = np.arange(0,8193,1024)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 4096):
        xticks_list = np.arange(0,4097,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
    elif (len(meas_data) == 2048):
        xticks_list = np.arange(0,2049,512)
        convxticks_list = data_converter(xticks_list, ecal_param)
        
    if convxticks_list[-1] <= 10000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/500) *500, 500)
    elif convxticks_list[-1] <= 25000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/1000) *1000, 1000)
    elif convxticks_list[-1] <= 50000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/2500) *2500, 2500)
    elif convxticks_list[-1] <= 100000:
        eV_xticks_list = np.arange(0, round(convxticks_list[-1]/5000) *5000, 5000)    
        
    meas_file_name = meas_file.split('.')[0]
    meas_data = meas_data[:-1]
    
    converted_x_axis = data_converter(np.linspace(0, meas_MCA-2, meas_MCA-1), ecal_param)
    
    fig, ax = plt.subplots(figsize=(10,4), dpi=DPI, nrows=2, sharex=True, height_ratios=[2,1])
    
    ax[0].plot(converted_x_axis, meas_data, color=color_schemes['c_five2'][0], label=f'Calibrated Spectrum \n ID: {meas_file_name}')
    
    for key, val in xraydb.xray_lines(elem_symbol).items():
        
        if both_axes == True:
            line_param_bothaxes = [val.intensity*np.max(meas_data),val.energy,1]
        line_param = [val.intensity,val.energy,1]
        # print(line_param[0:2])
        line_lin = np.linspace(val.energy-300,val.energy+300,1000)
        # xlin = np.linspace(0,int(E_max),int(E_max))
        
        if key[0] == 'K':
            col = color_schemes['c_rainbow'][1]
        elif (key[0] == 'L'):
            col = color_schemes['c_light'][1]
        elif (key[0] == 'M'):
            col = color_schemes['c_rainbow'][7]
        
        if both_axes == True:
            ax[0].plot(line_lin,gauss_func(line_param_bothaxes,line_lin), color=col, alpha=1, lw=1)
            # ax[0].plot(line_lin-1740,gauss_func(line_param_bothaxes,line_lin), color=col, alpha=0.5, lw=1) # Si KAlpha Escape Peak
        ax[1].plot(line_lin,gauss_func(line_param,line_lin), color=col, alpha=1, lw=1)
   
    ax[1].plot([],[],alpha=0,lw=1,label=f'Expected lines for {elem_symbol}')
    ax[1].plot([],[],color=color_schemes['c_rainbow'][1],lw=1,label='K-lines')
    ax[1].plot([],[],color=color_schemes['c_light'][1],lw=1,label='L-lines')
    ax[1].plot([],[],color=color_schemes['c_rainbow'][7],lw=1,label='M-lines')
    
    ax[0].set_xlim(0, convxticks_list[-1])
    ax[0].set_xticks(eV_xticks_list, [f'{int(i)}' for i in eV_xticks_list])
    
    if log_flag == True:
        ax[0].set_yscale('log')
        ax[0].set_ylim(1e0,np.max(meas_data)*1.2)
    else:
        ax[0].set_yscale('linear')
        ax[0].set_ylim(0,np.max(meas_data)*1.2)
    
    for i in range(2):
        ax[i].legend()
        ax[i].grid()
    
    plt.xlabel('Energy / eV', fontsize=10)
    ax[0].set_ylabel('Counts', fontsize=10)
    ax[1].set_ylabel('Intensity', fontsize=10)
 
    plt.tight_layout()
    
    if (both_axes == True):
        if (save_flag == True and log_flag == True):
            plt.savefig(f'./plots/calibrated/incl_expected_log/{meas_file_name}_expected_dualAx_log.png', transparent=False, dpi=DPI)
            plt.savefig(f'./plots/calibrated/incl_expected_log/{meas_file_name}_expected_dualAx_log.pdf', transparent=False, dpi=DPI)
        elif (save_flag == True and log_flag == False):
            plt.savefig(f'./plots/calibrated/incl_expected_lin/{meas_file_name}_expected_dualAx.png', transparent=False, dpi=DPI)
            plt.savefig(f'./plots/calibrated/incl_expected_lin/{meas_file_name}_expected_dualAx.pdf', transparent=False, dpi=DPI)
        else:
            plt.show()
    else:
        if (save_flag == True and log_flag == True):
            plt.savefig(f'./plots/calibrated/incl_expected_log/{meas_file_name}_expected_singleAx_log.png', transparent=False, dpi=DPI)
            plt.savefig(f'./plots/calibrated/incl_expected_log/{meas_file_name}_expected_singleAx_log.pdf', transparent=False, dpi=DPI)
        elif (save_flag == True and log_flag == False):
            plt.savefig(f'./plots/calibrated/incl_expected_lin/{meas_file_name}_expected_singleAx.png', transparent=False, dpi=DPI)
            plt.savefig(f'./plots/calibrated/incl_expected_lin/{meas_file_name}_expected_singleAx.pdf', transparent=False, dpi=DPI)
        else:
            plt.show()
    plt.close()
    return 0

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~#

def plot_data_multigauss_spectrum(meas_file:str, meas_MCA:int, meas_data:list, ecal_param:list, elem_symbol:str, save_flag:bool=False, log_flag:bool=False):
    
        
    return 11
        
d2026_03_25 = {
              '20260325-065918':'Cu', '20260325-070217':'Cu', '20260325-070325':'Cu',
              '20260325-080539':'Rb', '20260325-081037':'Rb', '20260325-081432':'Rb', '20260325-081734':'Rb',
              '20260325-091953':'Mo', '20260325-092112':'Mo',
              '20260325-102319':'Ag', '20260325-102428':'Ag',
              '20260325-120842':'Ba', '20260325-121129':'Ba', '20260325-121419':'Ba',
              '20260325-131924':'Tb', '20260325-132058':'Tb'}
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
if __name__ == "__main__":
    start_routine = time.monotonic_ns()
    
    args = parser.parse_args()
    
    ecal_data = read_ecal('./energy_cali/09_Apr_2026_eCali.json')
    file_list = read_full_directory(args.path)
    
    saving_plots_flag = args.save
    log_plots_flag = args.log
    
    for f in tqdm.tqdm(range(len(file_list))):
        meas_file, meas_data, meas_MCA, meas_EG, ecal_curr_param, ecal_curr_errors = measurement_parameters(file_list[f], ecal_data)
    
        if args.simple_mode == 'single':
            plot_basic_singleAxis_spectrum(meas_file, meas_EG, meas_MCA, meas_data, ecal_curr_param, saving_plots_flag, log_plots_flag)
        elif args.simple_mode == 'dual':
            plot_basic_dualAxis_spectrum(meas_file, meas_EG, meas_MCA, meas_data, ecal_curr_param, saving_plots_flag, log_plots_flag)
        
        if args.detailed_mode == '1':
            plot_data_expected_singleAxis_spectrum(meas_file, meas_MCA, meas_data, ecal_curr_param,elem_symbol=d2026_03_25[file_list[f].split('//')[2].split('.')[0]],both_axes=False, save_flag=saving_plots_flag, log_flag=log_plots_flag)
        elif args.detailed_mode == '2':
            plot_data_expected_singleAxis_spectrum(meas_file, meas_MCA, meas_data, ecal_curr_param,elem_symbol=d2026_03_25[file_list[f].split('//')[2].split('.')[0]],both_axes=True, save_flag=saving_plots_flag, log_flag=log_plots_flag)
    
    print('---------------------------------------')
    end_routine = time.monotonic_ns()
    elapsed_routine = (end_routine - start_routine) / 1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.1f} s)')
    print('---------------------------------------')