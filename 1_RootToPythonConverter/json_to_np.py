import json
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def json_to_np(json_file_path, json_field):
    print(os.path.exists(json_file_path))
    
    with open(json_file_path, "r") as j:
        data = json.load(j)
        
    # print(data["RBS;1"])
    # Zugriff auf das gewünschte Feld (z. B. "mydata")
    array = np.array(data[json_field])

    print(array)
    return array

def array_to_histo(data):
    plt.figure(figsize=(4,4), dpi=250)
    plt.hist(data[data <= 5], bins=2000, color='firebrick', edgecolor='firebrick')
    plt.xlabel('Ekin')
    plt.ylabel('Counts')
    plt.yscale('linear')
    plt.ylim(0,10000)
    plt.title("Histogramm von Ekin")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Verwendung: python json_to_np.py input.json")
    else:
        data = json_to_np(sys.argv[1], sys.argv[2])
        array_to_histo(data)