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
import xraydb
import plotly
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from lmfit import Model
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from getmac import get_mac_address as gma
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
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

plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r"\usepackage{mathtools}",
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def dual_print(*args, **kwargs):
    print(*args)
    with open('bins_matching_log.txt', 'a') as out:
        print(*args, **kwargs, file=out)
        
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

mac = gma()
print('---------------------------------------')
print(f"Current MAC address: {mac}")
print('---------------------------------------')
# mac_dict = {'f4:b5:20:5e:ba:f2': ['C://Users//schum//Documents//Filing Cabinet//1_RootFilesGeant4', 'C://Users//schum//Documents//Filing Cabinet//2_jsonFiles'], # Office
#             '14:5a:fc:4f:e8:35': ['D://root_files_temp_storage', 'D://json_files_temp_storage'], # Laptop
#             '0x1a7dda7115'  : ['B://IBA//root', 'B://IBA//json']} # Home PC
# '''
# Each mac adresse leads to a pair of paths, the first being the folder, 
# where the root files are, the second, where the json files are supposed to be stored
# '''
# root_path = mac_dict[mac][0]
# json_path = mac_dict[mac][1]

color_schemes = load_colors()

end_setup = time.process_time_ns()
elapsed_setup = (end_setup - start_setup)/1e6

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} ms)')
print('---------------------------------------')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
def DepDepth(rho, volt):
    return np.sqrt(rho*volt)

def sqrt_func(x,a,b,c,d):
    return a*np.sqrt(b*x + c) + d

sqrtModel = Model(sqrt_func)
print(f'parameter names: {sqrtModel.param_names}')
print(f'independent variables: {sqrtModel.independent_vars}')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#





# filename1 = ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-53-37.h5"
# filename2 = ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-49-18.h5"
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
#TODO: collect all h5 of one detector
files_52148 = [
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-06-17.h5",
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-10-49.h5",
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-17-50.h5",
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-20-14.h5",
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-25-30.h5",
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-35-07.h5"
]
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def h5_data_extraction(h5File):
    with h5py.File(h5File, "r") as f:
        iv = f["IV_data"][:]
        voltage = iv["voltage"]
        current = iv["current"]
    return voltage, current

def h5_data_compactor(h5FileList):
    measurement_dict = {}
    for file in h5FileList:
        print(file)
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
        

def h5_plotter(h5dict):
    fig, ax = plt.subplots(figsize=(4,3), dpi=300)
    c_index = 0
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        
        ax.plot(measurement['voltage'],measurement['current']*10**6,lw=0.75,ls=':',marker='x',ms=2.75, color=color_schemes['c_complementary'][c_index],label=k)
        c_index += 1

    #----------------- Detector Image includer -----------------#
    img = plt.imread("PIIPS_3D_f.png")
    annotation = TextArea(f"PIIPS Detector \n ID: {measurement['det_id']}", textprops=dict(color="black", fontsize=5, multialignment='center'))
    imagebox = OffsetImage(img, zoom=0.05)
    stacked = VPacker(children=[imagebox, annotation],
                 align="center",
                 pad=0,
                 sep=5)
    
    ab = AnnotationBbox(stacked, (50.,2.5), frameon=True)

    ax.add_artist(ab)
    #----------------- Detector Image includer -----------------#
    
    plt.xlabel(r'Bias Voltage / V')
    plt.ylabel(r'Leakage Current / $\mu$A')
    
    plt.xlim(0,500)
    # plt.ylim(1.8e-6,4e-6)
    plt.ylim(1.5,4)
    plt.grid(which='both')
    plt.legend(loc=2, fontsize=6)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    start_routine = time.process_time_ns()
    
    h5dict = h5_data_compactor(files_52148)
    h5_plotter(h5dict)
    
    print('---------------------------------------')
    end_routine = time.process_time_ns()
    elapsed_routine = (end_routine - start_routine)/1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.4f} ns)')
    print('---------------------------------------')