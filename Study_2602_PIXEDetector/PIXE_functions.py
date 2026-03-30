#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import csv
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
from scipy.odr import ODR, Model, RealData
from getmac import get_mac_address as gma
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
color_schemes = load_colors()
#----------------- Fitting Functions -----------------#
def sqrt_func(x,param):
    # return param[0]*np.sqrt(param[1]*x)
    return param[0]*np.sqrt(param[1]*x) + param[2]*x + param[3]
    # return param[0]*np.sqrt(param[1]*x - param[4]) + param[2]*x + param[3]

def lin_func(param, x):
    return param[0]*x + param[1]

def gauss_func(param, x):
    return param[0] * np.exp(-((x - param[1])**2)/(2*param[2]**2))

def double_gauss_func(param, x):
    return param[0] * np.exp(-((x - param[1])**2)/(2*param[2]**2)) + param[3] * np.exp(-((x - param[4])**2)/(2*param[5]**2))

def gauss_linear_func(param, x):
    G = gauss_func(param[0:3],x)
    L = lin_func(param[3:],x)
    return G + L

def energy_func(param, x):
    return param[0] + param[1]*x + param[2]*x**2
    
def exp_func(x, param):
    return param[0]*np.exp(param[1]*x) + param[2]

def cauchy_voigt_func(param, x):
    '''
    param[0]: sigma - std dev for Gaussian component \n
    param[1]: gamma - std dev for Lorentzian component \n
    '''
    
    return voigt_profile(x,param[0],param[1])








def evaluatorOLD(func, param_list:list, boundary_list:list, x:list, y:list, xerr:list, yerr:list):
    params = np.array(param_list)
    bounds_lower, bounds_upper = boundary_list[0], boundary_list[1]
    print(params)
    print(bounds_lower, bounds_upper)
        
    weight_x = 1/(np.array(xerr)**2)
    weight_y = 1/(np.array(yerr)**2)
    
    result = odr.odr_fit(func, x, y,
                         params, bounds=(bounds_lower, bounds_upper), 
                         weight_x=weight_x, weight_y=weight_y)
    
    print('STOP:', result.stopreason)
    print('PARAMS:', result.beta)
    
    return result.beta

def evaluator_scipy(func, beta0_list:list, x:list, y:list, xerr:list, yerr:list):
    
    data = RealData(x=x, y=y, sx=xerr, sy=yerr)
    model = Model(func)
    
    odr = ODR(data, model, beta0=beta0_list)
    output = odr.run()
    
    print('PARAMS:', output.beta)
    print('UNCERT:', output.sd_beta)
    
    return {'param':output.beta, 'errors':output.sd_beta}
#----------------- Fitting Functions -----------------#

def detector_pic(id):
    with open('RBS_detector_index.csv', newline='') as det_index:
        reader = csv.reader(det_index)
        index = list(reader)
    for d in range(len(index)):
        if index[d][0] == id:
            detector_type = index[d][1]
    if (detector_type == '1' or detector_type == 1):
        det_pic = 'SSB_3D_1.png'
    elif (detector_type == '2' or detector_type == 2):
        det_pic = 'SSB_3D_2.png'
    elif (detector_type == '3' or detector_type == 3):
        det_pic = 'SSB_3D_3_f.png'
    elif (detector_type == '4' or detector_type == 4):
        det_pic = 'SSB_3D_4_f.png'
    elif (detector_type == '5' or detector_type == 5):
        det_pic = 'PIIPS_3D_f.png'
    return det_pic

def peak_fitter(lines:list, peaks:list, peaks_err:list, init_values:list):
    
    lines_err = np.array([1]*len(lines))
    
    beta = evaluator_scipy(func=energy_func, beta0_list=init_values,
                           x=peaks, y=lines, xerr=peaks_err, yerr=lines_err)
    
    return beta


def read_json_formatted_file(filepath, encoding="utf-8"):
    """
    Reads a file whose contents are JSON-formatted, regardless of file extension.

    Parameters:
        filepath (str): Path to the file
        encoding (str): File encoding (default: utf-8)

    Returns:
        dict or list: Parsed JSON content

    Raises:
        ValueError: If the file content is not valid JSON
        OSError: If the file cannot be read
    """
    with open(filepath, "r", encoding=encoding) as f:
        content = f.read()

    try:
        return json.loads(content)
    except json.JSONDecodeError as e:
        raise ValueError(f"File content is not valid JSON: {e}") from e
    
def file_collector(measurement:str):
    file_collection = []
    top_level_path = f'.//{measurement}'
    for file in os.listdir(top_level_path):
        if (file[-5:] == '.vspc'):
            full_path = top_level_path + '//' + file
            file_collection.append(full_path)
    return file_collection

def pixe_single_spectrum_plot(filename:str):
    '''
    This function will produce a simple labeled plot of uncalibrated raw-data.
    '''
    DPI = 250
    data = read_json_formatted_file(filename)
    meas_name = filename.split('//')[2].split('.')[0]
    meas_folder = filename.split('//')[1]
    
    energyPerBin = data['Calibration']['BinSize_keV/Bin'] # keV/bin
    bin_data = data['RawData'][:-1] #remove that overflow bin at position 8191
    total_counts = np.sum(bin_data)
    total_counts_incl = np.sum(data['RawData'])
    # print(bin_data)
    
    bins = np.arange(0,len(bin_data),1)
    
    scatter_color = color_schemes['c_dark']
    
    fig, ax = plt.subplots(figsize=(6,3), dpi=DPI)
    # ax.set_facecolor(color_schemes['c_back'])
    # ax.plot(bins, bin_data, lw=0.75, color=scatter_color[0], zorder=2)
    ax.step(bins, bin_data, lw=0.75, color=scatter_color[3], zorder=2)
    
    if (len(bin_data) == 8191):
        ax.set_xlim(0,8193)
        ax.set_xticks(np.arange(0,8193,1024),np.arange(0,8193,1024))
    elif (len(bin_data) == 4095):
        ax.set_xlim(0,4095)
        ax.set_xticks(np.arange(0,4097,512),np.arange(0,4097,512))
    elif (len(bin_data) == 2047):
        ax.set_xlim(0,2047)
        ax.set_xticks(np.arange(0,2049,512),np.arange(0,2049,512))
        
        
    plt.xlabel('MCA channel')
    plt.ylabel('Counts')
    
    plt.grid(which="both")
    plt.tight_layout()

    #----------------- Information Box -----------------#
    #det_pic_file = detector_pic(measurement['det_id'])
    # img = plt.imread(det_pic_file)
    annotation = TextArea(f"X-ray measurement \n {meas_name} \n RAW DATA \n Total Counts: {total_counts}", textprops=dict(color="black", fontsize=5, multialignment='center'))
    # imagebox = OffsetImage(img, zoom=0.05)
    stacked = VPacker(children=[annotation],
                 align="center",
                 pad=0,
                 sep=5)
    
    ab = AnnotationBbox(offsetbox=stacked, xy=(0.9,0.85), xycoords='axes fraction', frameon=True)

    ax.add_artist(ab)
    #----------------- Information Box -----------------#

    
    plt.savefig(f'./plots/uncalibrated/{meas_folder}/{meas_name}.png', dpi=DPI)
    plt.savefig(f'./plots/uncalibrated/{meas_folder}/{meas_name}.pdf', dpi=DPI)
    
    plt.show()

def all_files_from_measSet(m_name:str):
    f_c = file_collector(m_name)
    for file in f_c:
        pixe_single_spectrum_plot(file)
    return 4