#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import h5py
import math
import tqdm
import xraydb
import plotly
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import odrpack as odr
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

def lam(m1, m2):
    return m1/m2

def sig(m1, m2, theta):
    '''
    INPUTS:\n
        m1, m2: amu \n
        theta: deg
    '''
    
    theta = theta * np.pi() / 180
    return 1 - lam(m1, m2)**2 * np.sin(theta)**2

def kap(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    theta = theta * np.pi() / 180
    return lam(m1, m2) * np.cos(theta)

def K_factor(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    return ((np.sqrt(sig(m1, m2, theta)) + kap(m1, m2, theta))/(1 + lam(m1, m2)))**2