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
from matplotlib.ticker import LinearLocator
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

def threeDimensional_k_factor_plot(m1:float, m2:list, theta:list):
    
    m1_m2 = []
    m2_alone = []
    for i in m2:
        m1_m2.append(m1/i.atomic_weight)
        m2_alone.append(i.atomic_weight)

    
    fig, ax = plt.subplots(figsize=(6,3), dpi=250, subplot_kw={"projection":"3d"})
    
    X = m1_m2
    Y = theta
    X,Y = np.meshgrid(X,Y)
    Z = rbs.kinematic_factor.K_factor(m1, m2_alone, Y)
    
    surface = ax.plot_surface(X,Y,Z, cmap='plasma', linewidth=0, antialiased=False)
    
    ax.set_zlim=(0,1)
    
    
    # Add a color bar which maps values to colors.
    fig.colorbar(surface, shrink=0.75, aspect=12)
    ax.set_box_aspect((16, 16, 8))
    ax.set_xlabel(r'Mass ration $\lambda$')
    ax.set_ylabel(r'Angle $\theta$')
    ax.set_zlabel(r'Kinematic factor $K$')
    
    plt.show()
    return 1

if __name__ == "__main__":
    m1 = md.element('He').atomic_weight
    m2 = [md.element(i) for i in range(2,83)]
    theta = np.arange(0,180,1)
    threeDimensional_k_factor_plot(m1, m2, theta)