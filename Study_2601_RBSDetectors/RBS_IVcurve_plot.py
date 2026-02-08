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

def sqrt_func(x,a,b,c):
    return a*np.sqrt(np.abs((x - b))) + c

       
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#


# filename1 = ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-53-37.h5"
# filename2 = ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-49-18.h5"
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
#TODO: collect all h5 of one detector
files_52148 = [
    ".//PIIPS//52148//IVcurve_PIIPS_52148___degC__2026-01-21_15-40-02.h5", #0-30
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2degC__2026-01-23_10-58-03.h5", #0-50
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-06-17.h5", #50-200
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-10-49.h5", #160-300
    ".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-25-30.h5" #300-460
]

files_57341 = [
    ".//PIIPS//57431//IVcurve_PIIPS_57431___dark_v2_noOscidegC__2026-01-26_15-49-19.h5",
    ".//PIIPS//57431//IVcurve_PIIPS_57431___dark_v2_noOscidegC__2026-01-23_13-27-00.h5"
]

files_57342 = [
    ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-49-18.h5",
    ".//PIIPS//57432//IVcurve_PIIPS_57432___dark_v2_noOscidegC__2026-01-23_13-53-37.h5"
]

files_33_268B = [
    ".//SSB//33-268B//IVcurve_SSB_33-268B___TC-021-300-300degC__2026-01-26_16-49-00.h5",
    ".//SSB//33-268B//IVcurve_SSB_33-268B___TC-021-300-300degC__2026-01-26_16-53-01.h5",
    ".//SSB//33-268B//IVcurve_SSB_33-268B___TC-021-300-300degC__2026-01-26_16-58-47.h5"
]

files_29_286 = [
    ".//SSB//29-286//IVcurve_SSB_29-286___MH-21-450-100degC__2026-01-26_16-26-04.h5",
    ".//SSB//29-286//IVcurve_SSB_29-286___MH-21-450-100degC__2026-01-26_16-22-37.h5"
]

files_29_286_single = [
    ".//SSB//29-286//IVcurve_SSB_29-286___MH-21-450-100degC__2026-01-26_16-26-04.h5"
]


# files_52148_short = [".//PIIPS//52148//IVcurve_PIIPS_52148___dark_v2_noOscidegC__2026-01-23_11-10-49.h5"]
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
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

def h5_measurement_combiner(h5dict):
    measure_total_voltage, measure_total_current = [],[]
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        for j in range(len(measurement['voltage'])):
            measure_total_voltage.append(float(measurement['voltage'][j]))
            measure_total_current.append(float(measurement['current'][j]))
    return np.array(measure_total_voltage), np.array(measure_total_current)

def h5_plotter(h5dict):
    DPI = 300
    plot_type = 'Zero-To-Breakdown'
    fig, ax = plt.subplots(figsize=(6,3), dpi=DPI)
    c_index = 0
    x_min, x_max, y_min, y_max = 1e6,0,1e6,0
    
    #----------------- Measurement Plotter -----------------#
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        scaling = 1e6
        
        if (x_min > np.min(measurement['voltage'])):
            x_min = np.min(measurement['voltage'])
            
        if (x_max < np.max(measurement['voltage'])):
            x_max = np.max(measurement['voltage'])
            
        if (y_min > np.min(measurement['current'])):
            y_min = np.min(measurement['current'])
        
        if (y_max < np.max(measurement['current'])):
            y_max = np.max(measurement['current'])
            
        # p_list =[0.5,np.min(measurement['voltage']),float(np.min(measurement['current']))*scaling]
        # print(p_list)
        # s_list = [np.min(measurement['voltage'])*0.95, np.max(measurement['voltage'])*1.1, 1000]
        
        ax.errorbar(x=measurement['voltage'],
                    y=measurement['current']*scaling,
                    yerr=0.05,
                    capsize=1,
                    capthick=0.4,
                    elinewidth=0.2,
                    ms=3,
                    fmt='x',
                    color=color_schemes['c_complementary'][c_index],
                    label=k)
        
        # result = evaluator(sqrtModel=sqrtModel, param_list=p_list, space=s_list, x=measurement['voltage'], y=scaled_measurement)
        
        # ax.plot(measurement['voltage'],
        #         result.best_fit,lw=0.75,
        #         ls='-',
        #         color=color_schemes['c_complementary'][c_index],
        #         alpha=0.9,
        #         zorder=2,
        #         label=k)
        
        c_index += 1
    #----------------- Measurement Plotter -----------------#
    
    #----------------- Detector Image includer -----------------#
    # img = plt.imread("PIIPS_3D_f.png")
    img = plt.imread("SSB_3D_3_f.png")
    annotation = TextArea(f"{measurement['det_type']} Detector \n ID: {measurement['det_id']}", textprops=dict(color="black", fontsize=5, multialignment='center'))
    imagebox = OffsetImage(img, zoom=0.05)
    stacked = VPacker(children=[imagebox, annotation],
                 align="center",
                 pad=0,
                 sep=5)
    
    ab = AnnotationBbox(stacked, (10.,0.8), frameon=True)

    ax.add_artist(ab)
    #----------------- Detector Image includer -----------------#
       
    #----------------- Plot Finalization -----------------#
    plt.xlabel(r'Bias Voltage / V')
    plt.ylabel(r'Leakage Current / $\mu$A')
    # plt.ylabel(r'Leakage Current / nA')
    
    if (x_min < 5 and x_max >= 100):
        x_min = -20
    
    plt.xlim(x_min*0.5, x_max*1.05)
    # plt.xlim(x_min*0.9, 375)
    # plt.ylim((y_min*scaling)*0.5, (y_max*scaling)*1.04)
    # plt.ylim((y_min*scaling)*0.5, (y_max*scaling)*1.35)
    plt.ylim(0.007,4)
    
    # plt.xlim(0,500)
    # plt.ylim(1.8e-6,4e-6)
    # plt.ylim(1.5,4)
    
    plt.yscale('log')
    
    plt.grid(which='both')
    plt.legend(loc='lower right', fontsize=6)
    plt.tight_layout()
    
    
    save_name_PIIPS = f'{measurement['det_type']}_' + plot_type + '_' + f'{measurement['det_id']}'
    plt.savefig('.//plots//' + save_name_PIIPS + '.pdf', dpi=DPI, transparent=False)
    plt.savefig('.//plots//' + save_name_PIIPS + '.png', dpi=DPI, transparent=False)
    
    plt.show()
    #----------------- Plot Finalization -----------------#
    return 0


if __name__ == "__main__":
    start_routine = time.process_time_ns()
    
    # h5dict = h5_data_compactor(files_52148)
    # h5dict = h5_data_compactor(files_57341)
    # h5dict = h5_data_compactor(files_57342)
    # h5dict = h5_data_compactor(files_33_268B)
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