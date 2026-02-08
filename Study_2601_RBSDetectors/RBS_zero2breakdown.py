#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.process_time_ns()
print('---------------------------------------')
print(time.strftime("RBS_zero2breakdown.py started: %a, %d %b %Y %H:%M:%S", time.localtime()))
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
color_set = color_schemes['c_rainbow']

end_setup = time.process_time_ns()
elapsed_setup = (end_setup - start_setup)/1e6

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} ms)')
print('---------------------------------------')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
def plot_zero2breakdown(h5dict):
    DPI = 300
    plot_type = 'Zero-To-Breakdown'
    
    fig, ax = plt.subplots(figsize=(6,3), dpi=DPI)
    c_index = 1
    x_min, x_max, y_min, y_max = 1e6,0,1e6,0
    
    #----------------- Measurement Plotter -----------------#
    for k in list(h5dict.keys()):
        measurement = h5dict[k]
        if (measurement['det_type'] == 'SSB'):
            scaling = 1e6
        else:
            scaling = 1e9
        
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
                    color=color_set[c_index],
                    label=k)
        
        # result = evaluator(sqrtModel=sqrtModel, param_list=p_list, space=s_list, x=measurement['voltage'], y=scaled_measurement)
        
        # ax.plot(measurement['voltage'],
        #         result.best_fit,lw=0.75,
        #         ls='-',
        #         color=color_schemes['c_complementary'][c_index],
        #         alpha=0.9,
        #         zorder=2,
        #         label=k)
        if (len(list(h5dict.keys())) == 2):
            c_index += 5
        elif (len(list(h5dict.keys())) == 3):
            c_index += 2
        elif (len(list(h5dict.keys())) > 3):
            c_index += 1
    #----------------- Measurement Plotter -----------------#
    
    save_name_PIIPS = f'{measurement['det_type']}_' + plot_type + '_' + f'{measurement['det_id']}'
    
    #----------------- Detector Image includer -----------------#
    det_pic_file = detector_pic(measurement['det_id'])
    img = plt.imread(det_pic_file)
    annotation = TextArea(f"{measurement['det_type']} Detector \n ID: {measurement['det_id']}", textprops=dict(color="black", fontsize=5, multialignment='center'))
    imagebox = OffsetImage(img, zoom=0.05)
    stacked = VPacker(children=[imagebox, annotation],
                 align="center",
                 pad=0,
                 sep=5)
    
    ab = AnnotationBbox(offsetbox=stacked, xy=(0.07,0.8), xycoords='axes fraction', frameon=True)

    ax.add_artist(ab)
    #----------------- Detector Image includer -----------------#
       
    #----------------- Plot Finalization -----------------#
    plt.xlabel(r'Bias Voltage / V')
    
    if (measurement['det_type'] == 'SSB'):
        plt.ylabel(r'Leakage Current / $\mu$A')
    else:
        plt.ylabel(r'Leakage Current / nA')
    
    if (x_min < 5 and x_max >= 100):
        x_min = -20
    
    plt.xlim(x_min*0.5, x_max*1.05)
    plt.ylim((y_min*scaling)*0.5, (y_max*scaling)*1.04)
    # plt.ylim((y_min*scaling)*0.5, (y_max*scaling)*1.35)
    
    if (len(sys.argv) > 3 and sys.argv[3] == 'log'):
        plt.yscale('log')
        plt.ylim(1e-1, (y_max*scaling)*1.04)
        save_name_PIIPS += '_log'
    else:
        pass
    
    plt.grid(which='both')
    plt.legend(loc='lower right', fontsize=6)
    plt.tight_layout()
    
    
    
    plt.savefig('.//plots//' + save_name_PIIPS + '.pdf', dpi=DPI, transparent=False)
    plt.savefig('.//plots//' + save_name_PIIPS + '.png', dpi=DPI, transparent=False)
    
    plt.show()
    #----------------- Plot Finalization -----------------#
    return 0


if __name__ == "__main__":
    start_routine = time.process_time_ns()
    
    fc = file_collector(dettype=sys.argv[1],id=sys.argv[2])
    h5dict = h5_data_compactor(fc)
    
    plot_zero2breakdown(h5dict)
    
    print('---------------------------------------')
    end_routine = time.process_time_ns()
    elapsed_routine = (end_routine - start_routine)/1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.4f} s)')
    print('---------------------------------------')
    
'''
prompt: py -3.13 RBS_zero2breakdown.py <det type> <det id> <'' or log>
'''