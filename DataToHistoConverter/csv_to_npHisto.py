import csv
import json
import sys
import os
import numpy as np
import pandas as pd

def csv_to_npHisto(file):
    # print(os.path.exists(file))
    df = pd.read_csv(file, delimiter=';',header=7)
    histo_values = []
    for i in range(len(df)):
        sub = np.array([float(df.at[i, 'Energy [eV]']*10**(-6))]*int(df.at[i, 'Counts']))
        # print(sub)
        histo_values += list(sub)
    # print(histo_values)
    return histo_values

if __name__ == "__main__":
    csv_to_npHisto(sys.argv[1])