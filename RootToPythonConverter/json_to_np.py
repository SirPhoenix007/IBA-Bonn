import json
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def json_to_np(json_file_path, json_field=None):
    print(os.path.exists(json_file_path))
    if (json_field != None):
        with open(json_file_path, "r") as j:
            data = json.load(j)
            
        # print(data["RBS;1"])
        # Zugriff auf das gewünschte Feld (z. B. "mydata")
        if (json_field != 'all'):
            array = np.array(data[json_field])

            print(array)
            return array
        else:
            array_set = {}
            i = 1
            for k in data.keys():
                array = np.array(data[k])
                array_set[k] = array
                print(f'Array for key {k} created. {i}/{len(data.keys())}')
                i += 1
            print(f'WARNING: The data is now structured in a dictionary, with the keys {data.keys()}')
            return array_set
    else:
        print('ERROR: No json field specified!')
        return 41
    

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
        # array_to_histo(data)