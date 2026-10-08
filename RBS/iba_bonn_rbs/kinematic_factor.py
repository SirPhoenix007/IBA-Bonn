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
        theta: rad
    '''
    return 1 - lam(m1, m2)**2 * np.sin(theta)**2

def kap(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: rad
    '''
    return lam(m1, m2) * np.cos(theta)

def K_factor(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    theta = theta * np.pi / 180
    
    return ((np.sqrt(sig(m1, m2, theta)) + kap(m1, m2, theta))/(1 + lam(m1, m2)))**2

def dKdt(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    theta = theta * np.pi / 180
    
    return K_factor(m1, m2, theta)*(-2*lam(m1, m2)*np.sin(theta))/(np.sqrt(sig(m1, m2, theta)))

def dKdM1(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    theta = theta * np.pi / 180
    
    A = (np.cos(theta)*(2*sig(m1,m2,theta)-1)*np.sqrt(sig(m1,m2,theta))**(-1) + lam(m1,m2)*np.cos(2*theta)*np.sqrt(sig(m1, m2, theta)))/((1+lam(m1,m2))**2)
    B = K_factor(m1, m2, theta)/(1+lam(m1, m2))
    
    return 2/m2 * (A - B)

def dKdM2(m1, m2, theta):
    '''
    INPUTS: \n
        m1, m2: amu \n
        theta: deg
    '''
    theta = theta * np.pi / 180
    
    A = (lam(m1, m2)**2*np.cos(2*theta) - (kap(m1, m2, theta))/(np.sqrt(sig(m1, m2, theta))))/((1+lam(m1,m2))**2)
    B = lam(m1, m2)*K_factor(m1, m2, theta)/(1+lam(m1, m2))
    
    return 2/m2 * (A - B)