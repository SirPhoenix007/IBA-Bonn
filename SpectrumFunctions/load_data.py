import numpy as np
import json
import csv
import pandas as pd

from RootToPythonConverter.json_to_np import *
from DataToHistoConverter.csv_to_npHisto import *

def Load_Data(file_name, opt = None):
    suffix = file_name.split('.')[-1]
    if (suffix == 'csv'):
        if (opt == None):
            data_df = pd.read_csv(file_name, delimiter=';',header=7)
            return data_df
        elif (opt == 'histo'):
            data_histo = csv_to_npHisto(file_name)
            return data_histo
    if (suffix == 'json'):
        data = json_to_np(file_name , json_field='all')
        return data