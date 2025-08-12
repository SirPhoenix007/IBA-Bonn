import json
import sys
import os
import uuid
import numpy as np

from root_to_json import *
from json_to_np import *

mac = hex(uuid.getnode())
print(f"MAC address: {mac}")

mac_dict = {'0xb6ab0b4445f9': ['C://Users//schum//Documents//Filing Cabinet//1_RootFilesGeant4', 'C://Users//schum//Documents//Filing Cabinet//2_jsonFiles'], # Office
            '0x145afc4fe836': ['C://Users//schum//Documents//root_files_temp_storage', 'C://Users//schum//Documents//json_files_temp_storage'], # Laptop
            '0x1a7dda7115': ['B://IBA//root', 'B://IBA//json']} # Home PC
# Each mac adresse leads to a pair of paths, the first being the folder, 
# where the root files are, the second, where the json files are supposed to be stored
root_path = mac_dict[mac][0]
json_path = mac_dict[mac][1]

print(root_path, json_path)
files_r, files_j = [],[]
convert_size = 0

print('ROOT-Files to be converted to JSON:') 
for rootfile in os.listdir(root_path):
    found = False
    root_file_path = root_path + '//' + rootfile + '//Root_Files//'
    
    file_name = next(f for f in os.listdir(root_file_path) if os.path.isfile(os.path.join(root_file_path, f)))
    
    full_root_file_path = os.path.join(root_file_path, file_name)
    full_json_file_path = json_path + '//' + file_name[:-5]
    
    for f in os.listdir(json_path):
        if f.startswith(file_name[:-5]):
            found = True
    # print(full_json_file_path)
    if (found == False):
        files_r.append(full_root_file_path)
        files_j.append(full_json_file_path)
        convert_size += os.path.getsize(full_root_file_path)/(10**6)
        print(full_root_file_path)

if (len(files_r) != 0):
    print('--------------------------------------------------------')
    print(f'The number of ROOT-files needed to be converted is: {len(files_r)}')
    print('--------------------------------------------------------')
    print(f'Total size of these ROOT-files is {round(convert_size,2)} MB.')
    pass
else:
    print('All current ROOT-files are already converted.')

# root_to_json(files_r[7], files_j[7])
for f in range(len(files_r)):
    root_to_json(files_r[f], files_j[f])


    