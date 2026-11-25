
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
from colors import load_colors
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
cs = color_schemes['c_complementary']

plt.rcParams["axes.grid"] = True


def lin_func(x, a, b):
    return a * x + b

lines = np.array([277,524.9,3311.1,3313.8])
peak_error = 10

holm_bean_area_peaks    = [280,520,3320,3320]
holm_bean_peaks         = [280,530,3320,3320]
holm_ground_peaks       = [270,520,3330,3330]
rcvt_small_peaks        = [270,520,3320,3320]
rcvt_large_peaks        = [270,520,3310,3310]
rcvt_ground_peaks       = [270,520,3320,3320]

fig, ax = plt.subplots(figsize=(4,4), dpi=300, nrows=2, ncols=3, sharex=True, sharey=True)


popt, pcov = curve_fit(lin_func, lines, holm_bean_peaks)
ax[0,0].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[0,0].scatter(lines, holm_bean_peaks, color=cs[0], s=8)

popt, pcov = curve_fit(lin_func, lines, holm_bean_area_peaks)
ax[0,1].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[0,1].scatter(lines, holm_bean_area_peaks, color=cs[0], s=8)

popt, pcov = curve_fit(lin_func, lines, holm_ground_peaks)
ax[0,2].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[0,2].scatter(lines, holm_ground_peaks, color=cs[0], s=8)

popt, pcov = curve_fit(lin_func, lines, rcvt_small_peaks)
ax[1,0].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[1,0].scatter(lines, rcvt_small_peaks, color=cs[0], s=8)

popt, pcov = curve_fit(lin_func, lines, rcvt_large_peaks)
ax[1,1].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[1,1].scatter(lines, rcvt_large_peaks, color=cs[0], s=8)

popt, pcov = curve_fit(lin_func, lines, rcvt_ground_peaks)
ax[1,2].plot(lines, lin_func(lines, *popt), color=cs[1], lw=0.75, label=f'f(x) = {round(popt[1],3)} + {round(popt[0],3)} * x')
ax[1,2].scatter(lines, rcvt_ground_peaks, color=cs[0], s=8)

for x in ax.flat:
    x.legend(loc='best', fontsize=5)

fig.supxlabel('Database line energies in eV')
fig.supylabel('Measured line energies in eV')
plt.show()