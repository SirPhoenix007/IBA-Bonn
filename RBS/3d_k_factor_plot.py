#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import tqdm
import xraydb
import argparse
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
import odrpack as odr
import seaborn as sb
import mendeleev as md
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
from getmac import get_mac_address as gma
from itertools import chain
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import iba_bonn_rbs as rbs
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


plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def threeDimensional_k_factor_plot(m1:float, m2:list, theta:list, colors:list, linestyles:list):
    
    m2_atw = []
    m2_name = []
    for i in m2:
        m2_atw.append(i.atomic_weight)
        m2_name.append(i.name)
    
    plt.figure(figsize=(5,3), dpi=160)
    
    for m in range(len(m2)):
            k_list = []
            for t in theta:
                k_list.append(rbs.kinematic_factor.K_factor(m1, m2_atw[m], t))
                
            plt.plot(theta, k_list, color=colors[m], ls=linestyles[m], lw=1, label=f'Target = {m2_name[m]}')
        
    plt.xlabel(r'Angle $\theta$')
    plt.ylabel(r'Kinematic factor $K$')
    plt.grid()
    plt.legend()
    plt.show()
    return 1

if __name__ == "__main__":
    m1 = md.element('He').atomic_weight
    m2 = [md.element('O'),md.element('K'), md.element('Fe'), md.element('Ag'), md.element('Nd'), md.element('U')]
    theta = np.arange(0,180,1)
    colors = rbs.colors.load_colors()['c_dark']
    linestyles = ['-','-','-','-','-','-']
    threeDimensional_k_factor_plot(m1, m2, theta, colors, linestyles)