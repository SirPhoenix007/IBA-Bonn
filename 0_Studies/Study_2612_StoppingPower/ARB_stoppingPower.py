# Ashley, Ritchie & Brandt Stopping Power
#by H. Schumacher

lazy import argparse
lazy import numpy as np
lazy import matplotlib.pyplot as plt
lazy import time
lazy import tqdm
lazy from functools import wraps
lazy from scipy import integrate
lazy import scipy.special as scsp

'''
SOURCES:
Ashley_1972.pdf -> (1)
Fernandez_Varea_2026.pdf -> (2)
Sannikov_1958.pdf -> (3)
Jensen_1932.pdf -> (4)
Lenz_1932.pdf -> (5)
'''

#F-Function (1.14) / (2.13/14)
#Xi-Function
#x & w
#eta, chi
#L(x) function (1.13)

V_ZERO = 2.19*1e6

def x_parameter(Z1:int, Z2:int, E1:float):
    '''
    Index 1: Projectile\n
    Z1, E1 in MeV\n
    \n
    Index 2: Target\n
    Z2\n
    \n
    Conversion from (1)\n
    v_1^2/v_0^2 = 40.2*E1 in MeV/amu
    '''
    x = float((40.2*E1/Z1)/Z2)
    v1 = float(np.sqrt(40.2*E1/Z1 * V_ZERO**2))
    
    return x, v1


if __name__ == "__main__":
    print(x_parameter(2,13,1))