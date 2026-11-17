# RBS cross section calculation
# by Henry Schumacher

'''
RBS in CMS system
RBS Mayer -> arcsin(theta1/theta2)?!
non-RBS corrections
'''




lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import xraydb as xy
lazy import time
lazy import tqdm
lazy from functools import wraps
lazy from scipy import integrate
lazy from scipy.special import expn
# lazy import medeleev as md #TODO to be added in future, when medeleev is functional under 3.15

def rutherford_cross_section(E1, Z1, Z2, angle, verbose=False):
    '''
    LABORATORY FRAME
    E1: MeV
    Z1, Z2: int
    angle: degree
    '''
    if verbose:
        print('='*40)
        print(f'{xy.atomic_symbol(Z1)} on {xy.atomic_symbol(Z2)}')
        print('-'*40)
        print(f'Energy = {E1} MeV under \u03B8 =  {angle}°')
        print('='*40)
    
    
    M1 = xy.atomic_mass(Z1)
    M2 = xy.atomic_mass(Z2)
    
    theta = np.deg2rad(angle)
    lamb = M1 / M2

    ec = 1.44                   # MeV*fm
    
    prefactor = ((Z1*Z2*ec)/(4*E1))**2
    ruth_angle = 4 / (np.sin(theta)**4)
    numerator = (np.sqrt(1 - (lamb*np.sin(theta))**2) + np.cos(theta))**2
    denominator = np.sqrt(1 - (lamb*np.sin(theta))**2)
    
    result = prefactor * ruth_angle * numerator / denominator
    conversion = 1e-2
    return result*conversion

if __name__ == "__main__":
    # print(f'RESULT: {round(rutherford_cross_section(E1=5, Z1=2, Z2=79, angle=60)*10,2)} mb/sr ')
    print(rutherford_cross_section(5, 2, 79, 60))