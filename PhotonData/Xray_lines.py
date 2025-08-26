import xraydb
import json
import numpy as np
import os
import sys


def xray_lines_extractor():
    xray_lines = {}
    for Z in range(1, 101):  # Elements from Hydrogen (Z=1) to Fermium (Z=100)
        element = xraydb.atomic_symbol(Z)
        # Get all known X-ray lines for the element
        try:
            lines = xraydb.xray_lines(element)
        except Exception:
            continue
            
        if (lines != {}):
            xray_lines[element] = lines
    
    with open('xraylines.json', "w") as json_file:
        json.dump(xray_lines, json_file, indent=1)
        
        
        

if __name__ == "__main__":
    xray_lines_extractor()