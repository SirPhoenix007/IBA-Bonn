#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.process_time_ns()
print('---------------------------------------')
print(time.strftime("analysis.py started: %a, %d %b %Y %H:%M:%S", time.localtime()))
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import xraydb
import plotly
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import matplotlib.pyplot as plt
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

def dual_print(*args, **kwargs):
    print(*args)
    with open('bins_matching_log.txt', 'a') as out:
        print(*args, **kwargs, file=out)
        
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
#Functions Access
from DataToHistoConverter.csv_to_npHisto import *
from RootToPythonConverter.json_to_np import *
from RootToPythonConverter.colors import load_colors
import SpectrumFunctions as sf
# from SpectrumFunctions.spectra_class1 import *
# from SpectrumFunctions.load_data import *
# from SpectrumFunctions.spectra_class1 import *
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
#Photon Line Information
gamma_photon_line_energies_path = 'PhotonData//ENDSF_gamma_energy.json' #short dict to find lines via their energy
gamma_photon_line_fulldata_path = 'PhotonData//ENDSF_gamma_full.json' #long dict with same keys as above but with full information
xray_photon_lines_path          = 'PhotonData//xraylines.json'
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

mac = gma()
print('---------------------------------------')
print(f"Current MAC address: {mac}")
print('---------------------------------------')
mac_dict = {'f4:b5:20:5e:ba:f2': ['C://Users//schum//Documents//Filing Cabinet//1_RootFilesGeant4', 'C://Users//schum//Documents//Filing Cabinet//2_jsonFiles'], # Office
            '14:5a:fc:4f:e8:35': ['D://root_files_temp_storage', 'D://json_files_temp_storage'], # Laptop
            '0x1a7dda7115'  : ['B://IBA//root', 'B://IBA//json']} # Home PC
'''
Each mac adresse leads to a pair of paths, the first being the folder, 
where the root files are, the second, where the json files are supposed to be stored
'''
root_path = mac_dict[mac][0]
json_path = mac_dict[mac][1]

color_schemes = load_colors()

end_setup = time.process_time_ns()
elapsed_setup = (end_setup - start_setup)/1e6

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} ms)')
print('---------------------------------------')
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
def load_parameter(parameter_file):
    parameterFrame = pd.read_csv(parameter_file, header=0, delimiter=';')
    parameterDict = parameterFrame.to_dict(orient='records')
    print(f'LOADED PARAMETERS FROM {parameter_file}')
    print('---------------------------------------')
    print(parameterFrame)
    
    files, names = sf.file_merger.parameter_to_list(parameterDict)
    
    #SPECTRUM CLASS 1
    pmt = {}
    pmt['col_s']    = parameterDict[9]['VALUE'] # color scheme
    pmt['files']    = files # files as list
    pmt['names']    = names # spec. names as list
    pmt['plt_tlt']  = parameterDict[10]['VALUE'] # plot title
    
    sf.spectra_class1.spectrum_cl1(pmt)
    
    return 'done'

if __name__ == "__main__":
    start_routine = time.perf_counter_ns()
    load_parameter(sys.argv[1])
    
    print('---------------------------------------')
    end_routine = time.perf_counter_ns()
    elapsed_routine = (end_routine - start_routine)/1e9
    print(f'INFO: ROUTINE COMPLETE ({elapsed_routine:.4f} ns)')
    print('---------------------------------------')