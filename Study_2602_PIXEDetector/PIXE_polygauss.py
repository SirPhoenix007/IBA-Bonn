import numpy as np
import matplotlib.pyplot as plt

def gauss(center, height, width, x):
    return height * np.exp(-((x - center)**2)/(2*width**2))

def multi_gauss(param:list, x):
    total_gauss = np.zeros(len(x))
    for i in range(0,len(param),3):
        center = param[i]
        height = param[i+1]
        width  = param[i+2]
        total_gauss += gauss(center, height, width, x)
    return total_gauss