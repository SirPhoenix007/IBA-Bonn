
#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import time
import xraydb
import plotly
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from getmac import get_mac_address as gma
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
def dual_print(output, txt):
    print(txt)
    with open(output, 'a') as out:
        print(txt, file=out)
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import SpectrumFunctions as sf
from RootToPythonConverter.json_to_np import *
from DataToHistoConverter.csv_to_npHisto import *
from RootToPythonConverter.colors import load_colors
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#


mac = gma()


mac_dict = {'f4:b5:20:5e:ba:f2': ['C://Users//schum//Documents//Filing Cabinet//1_RootFilesGeant4', 'C://Users//schum//Documents//Filing Cabinet//2_jsonFiles'], # Office
            '14:5a:fc:4f:e8:35': ['D://root_files_temp_storage', 'D://json_files_temp_storage'], # Laptop
            '0x1a7dda7115'  : ['B://IBA//root', 'B://IBA//json']} # Home PC
# Each mac adresse leads to a pair of paths, the first being the folder, 
# where the root files are, the second, where the json files are supposed to be stored
# root_path = mac_dict[mac][0]
json_path = mac_dict[mac][1]

color_schemes = load_colors()

# time.strftime("Kernel started: %a, %d %b %Y %H:%M:%S", time.localtime())

def spectrum_cl1(parameters_cl1:dict):
    
    '''
    INPUTS:\n
    parameter_cl1: dictionary
    '''
    timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    
    #PARSING OF PARAMETERS FOR PLOT
    color_scheme    = parameters_cl1['col_s']
    files           = parameters_cl1['files']
    dictnames       = parameters_cl1['names']
    plot_title      = parameters_cl1['plt_tlt']
    save_fig        = parameters_cl1['savefig']
    fig_name        = parameters_cl1['fig_name']
    
    DPI = 300
    
    colors = color_schemes[color_scheme]
    plt.figure(figsize=(7,4), dpi=DPI)
    
    print('---------------------------------------')
    print('RUNNING NOW: Spectrum Class 1')
    print('---------------------------------------')
    
    for file in range(len(files)):
        print(files[file])
        df = sf.load_data.Load_Data(file_name=files[file])
        x = df['Energy [eV]']
        y = df['Counts']
        plt.plot(x, y,
                 lw=0.75, ls='-',
                 color=colors[file], label= dictnames[file],
                 zorder=2)
    plt.grid()
    plt.legend(fontsize=6)
    
    plt.xscale('linear')
    plt.xlim(0,4000)
    plt.xlabel('Energy in eV')
    
    plt.yscale('log')
    plt.ylim(100,2*10**5)
    plt.ylabel('Counts')
    
    plt.title(plot_title)
    plt.tight_layout()
    if (save_fig == 'True'):
        temp_file_path = files[file].split('//')[:-1]
        file_path = ''
        for i in range(len(temp_file_path)):
            file_path += f'{temp_file_path[i]}//'
        save_name = 'plot_' + fig_name + '_' + timestamp
        plt.savefig(file_path + save_name + '.pdf', dpi=DPI, transparent=False)
        plt.savefig(file_path + save_name + '.png', dpi=DPI, transparent=False)
    
    plt.show()
    