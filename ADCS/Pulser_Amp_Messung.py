#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
start_setup = time.time()
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
    with open('xxx.txt', 'a') as out:
        print(*args, **kwargs, file=out)
        
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

color_schemes = load_colors()

end_setup = time.time()
elapsed_setup = (end_setup - start_setup)

print(f'INFO: SETUP COMPLETE ({elapsed_setup:.2f} s)')
print('---------------------------------------')

data = np.loadtxt('Pulser_Amp_010726.txt', delimiter=',')
# print(data)

plt.figure(figsize=(6,4), dpi=300)

plt.plot(data[:,0],data[:,1]*10**(-3), color=color_schemes['c_complementary'][0], zorder=2, label='Pulser Signal')
plt.plot(data[:,0],data[:,2], color=color_schemes['c_complementary'][2], zorder=2, label='Amplifier Signal')
plt.plot(data[:,0],data[:,3], color=color_schemes['c_complementary'][3], zorder=2, label='Linear Gate Stretcher Signal')

plt.xlabel('Pulser Potentiometer Setting')
plt.ylabel('Pulser Signal Height / V')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,4), dpi=300)

plt.plot(data[:,0],data[:,4], color=color_schemes['c_complementary'][3], zorder=2, label=r'Time Amp. $->$ LGS')
plt.plot(data[:,0],data[:,5], color=color_schemes['c_complementary'][4], zorder=2, label=r'Time Pulse $->$ Amp.')


plt.xlabel('Pulser Potentiometer Setting')
plt.ylabel(r'Delay / $\mu$s')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()