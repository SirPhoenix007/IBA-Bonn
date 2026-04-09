
#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.process_time_ns()
print('---------------------------------------')
print(time.strftime("PIXE_plotting.ipynb started: %a, %d %b %Y %H:%M:%S", time.localtime()))
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
from scipy.special import voigt_profile
from getmac import get_mac_address as gma
from itertools import chain
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
from PIXE_functions import *
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
    print(f'INFO: FILE {file} READ (MCA: {curr_MCA}, EG: {curr_EG})')
    
    for i in range(1,6):
        if curr_MCA == ecal_data[str(i)]['MCA'] and curr_EG == ecal_data[str(i)]['EG']:
            ecal_MCA = ecal_data[str(i)]['MCA']
            ecal_EG = ecal_data[str(i)]['EG']
            ecal_param = ecal_data[str(i)]['param']
            ecal_paramUn = ecal_data[str(i)]['paramUn']
            print(f'INFO: ECAL for FILE: DC-offset {ecal_param[0]:.4f} | Linear {ecal_param[1]:.4f} | Quadratic {ecal_param[2]*1e6:.4f} x 10^-6')
            continue
    return curr_raw_data, curr_MCA, curr_EG, ecal_param, ecal_paramUn

def data_converter(measurement_data:list, ecal_param:list):
    converted_data = []
    for i in range(len(measurement_data)):
        converted_value = ecal_param[0] + ecal_param[1]*measurement_data[i] + ecal_param[2]*measurement_data[i]**2
        converted_data.append(converted_value)
    return converted_data

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~#

def plot_basic_spectrum(meas_EG:int, meas_MCA:int, meas_data:list, ecal_param:list):
    converted_data = data_converter(meas_data, ecal_param)
    converted_x_axis = data_converter(np.linspace(0, meas_MCA-1, meas_MCA), ecal_param)
    plt.figure(figsize=(6,4))
    plt.plot(converted_x_axis, converted_data, color=color_schemes['c_five2'][0], label='Calibrated Spectrum')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Counts')
    plt.title(f'Calibrated Spectrum (EG: {meas_EG} eV, MCA: {meas_MCA})')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
if __name__ == "__main__":
    start_routine = time.monotonic_ns()
    ecal_data = read_ecal('./energy_cali/09_Apr_2026_eCali.json')
    file_list = read_full_directory('collected_data//2026_03_25')
    meas_data, meas_MCA, meas_EG, ecal_curr_param, ecal_curr_errors = measurement_parameters(file_list[-1], ecal_data)
    
    plot_basic_spectrum(meas_EG, meas_MCA, meas_data, ecal_curr_param)
    
    print('---------------------------------------')
    end_routine = time.monotonic_ns()
    elapsed_routine = (end_routine - start_routine) / 1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.1f} s)')
    print('---------------------------------------')