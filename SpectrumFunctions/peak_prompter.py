
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
def xray_line_searcher(energy_low, energy_high):
    '''
    energy_low/high: in eV
    '''
    matching = []
    xray_line_path = './/xraylines.json'
    
    if (xraydb.xray_lines('Fe') != {}):
        for Z in range(1, 101):  # Elements from Hydrogen (Z=1) to Fermium (Z=100)
            element = xraydb.atomic_symbol(Z)
            # Get all known X-ray lines for the element
            try:
                lines = xraydb.xray_lines(element)
            except Exception:
                continue  # Skip elements that don't have data
            
            # print(lines)
            for line_id, line in lines.items():
                line_energy = line.energy  # Energy in keV
                if (energy_low <= line_energy and line_energy <= energy_high):
                    matching.append({
                        'element': element,
                        'Z': Z,
                        'transition': line_id,
                        'energy_eV': line_energy,
                        'intensity': line.intensity
                    })
    else: #the data extraction has to be fixed.
        with open(xray_line_path, "r") as j:
            xraydata = json.load(j)
        for Z in range(1,101):
            element = xraydb.atomic_symbol(Z)
            try:
                lines = xraydata[element]
            except Exception:
                continue
            
            for line_id, line in lines.items():
                line_energy = line.energy  # Energy in keV
                if (energy_low <= line_energy and line_energy <= energy_high):
                    matching.append({
                        'element': element,
                        'Z': Z,
                        'transition': line_id,
                        'energy_eV': line_energy,
                        'intensity': line.intensity
                    })

    return matching


def peak_text_prompter(peaks, energy_centers, energy_width, output_file_name, info):
    '''
    Energies from this function are in keV! \n
    '''
    energy_bins = {}
    ofn = output_file_name
    enw = energy_width
    enc = energy_centers
    
    dual_print(ofn,'--'.ljust(80,'-'))
    dual_print(ofn, info)
    dual_print(ofn,'--'.ljust(80,'-'))
    dual_print(ofn,'Bin'.ljust(8) + '| ' + 'E_range (eV)'.ljust(18) + '| ' + 'El.'.ljust(4)  + '| ' + 'Line'.ljust(6) + '| ' + 'Energy (eV)'.ljust(12) + '| ' + 'Int.'.ljust(6))
    
    for p in peaks:
        energy_bins[p] = [enc[p], enc[p] - enw, enc[p] + enw]
        
        E_low = round(energy_bins[p][1],4)
        E_high = round(energy_bins[p][2],4)
        matching = xray_line_searcher(E_low, E_high)
        # print(matching[0])
        # print(energy_bins[p])
        dual_print(ofn,f'Bin {p}'.ljust(8) + '| ' + f'{round(energy_bins[p][1],1)} - {round(energy_bins[p][2],1)}'.ljust(18) + '| ' + '-'.ljust(33,'-'))
        
        for line in matching:
            dual_print(ofn,''.ljust(28) + '| ' + str(line['element']).ljust(4) + '| ' + str(line['transition']).ljust(6) + '| ' + str(line['energy_eV']).ljust(12) + '| ' + str(round(line['intensity'],3)).ljust(6))
    
    # print(energy_bins)
    return energy_bins