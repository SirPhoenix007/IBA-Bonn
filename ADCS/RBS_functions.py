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
from getmac import get_mac_address as gma
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

#----------------- Fitting Functions -----------------#
def sqrt_func(x,param):
    # return param[0]*np.sqrt(param[1]*x)
    return param[0]*np.sqrt(param[1]*x) + param[2]*x + param[3]
    # return param[0]*np.sqrt(param[1]*x - param[4]) + param[2]*x + param[3]
    
def exp_func(x,param):
    return param[0]*np.exp(param[1]*x) + param[2]

def evaluator(func, param_list:list, boundary_list:list, x:list, y:list, xerr:list, yerr:list):
    params = np.array(param_list)
    bounds_lower, bounds_upper = boundary_list[0], boundary_list[1]
    print(params)
    print(bounds_lower, bounds_upper)
        
    weight_x = 1/(np.array(xerr)**2)
    weight_y = 1/(np.array(yerr)**2)
    
    result = odr.odr_fit(func, x, y,
                         params, bounds=(bounds_lower, bounds_upper), 
                         weight_x=weight_x, weight_y=weight_y)
    
    print(result.stopreason)
    print(result.beta)
    
    return result.beta
#----------------- Fitting Functions -----------------#

#----------------- h5-file Functions -----------------#
def file_collector(dettype:str,id:str):
    file_collection = []
    top_level_path = f'.//{dettype}//{id}'
    for file in os.listdir(top_level_path):
        if (file[-3:] == '.h5'):
            full_path = top_level_path + '//' + file
            file_collection.append(full_path)
    return file_collection


def h5_data_extraction(h5File):
    with h5py.File(h5File, "r") as f:
        iv = f["IV_data"][:]
        voltage = iv["voltage"]
        current = iv["current"]
    return voltage, current

def h5_data_compactor(h5FileList):
    measurement_dict = {}
    for file in h5FileList:
        # print(file)
        parts = file.split('//')
        det_type = parts[1]
        det_id = parts[2]
        det_type2 = parts[3].split('__')[1][1:-4]
        date = parts[3].split('__')[2][:-3]
        v,c = h5_data_extraction(file)
        c_min, c_max = np.min(c), np.max(c)
        measurement_dict[date] = {'det_type':det_type, 'det_type2':det_type2, 'det_id':det_id, 'voltage':v, 'current':c, 'c_bounds': (c_min,c_max)}
        # print(measurement_dict)
    return measurement_dict
        
def h5_measurement_combiner(h5dict):
    measure_total_voltage, measure_total_current = [],[]
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        for j in range(len(measurement['voltage'])):
            measure_total_voltage.append(float(measurement['voltage'][j]))
            measure_total_current.append(float(measurement['current'][j]))
    return np.array(measure_total_voltage), np.array(measure_total_current)
#----------------- h5-file Functions -----------------#

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