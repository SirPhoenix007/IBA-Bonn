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
from scipy.odr import ODR, Model, RealData
from getmac import get_mac_address as gma
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

#----------------- Fitting Functions -----------------#
def sqrt_func(x,param):
    # return param[0]*np.sqrt(param[1]*x)
    return param[0]*np.sqrt(param[1]*x) + param[2]*x + param[3]
    # return param[0]*np.sqrt(param[1]*x - param[4]) + param[2]*x + param[3]

def gauss_func(param, x):
    return param[0] * np.exp(-((x - param[1])**2)/(2*param[2]**2))
    
def exp_func(x, param):
    return param[0]*np.exp(param[1]*x) + param[2]

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