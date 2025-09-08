import numpy as np
import json
import csv
import pandas as pd

def parameter_to_list(param_dict):
    files = []
    names = []
    i = 1
    while i != 9:
        if (type(param_dict[i]['VALUE']) == str):
            files.append(param_dict[i]['VALUE'])
            names.append(param_dict[i]['COMMENT'])
        else:
            break
        i += 1
    return files, names
        