import pandas as pd
import csv
import numpy as np
import json
import os
import sys


#ENERGIES ARE IN KEV

def file_analyser():
    
    DIR = './/gamma_data'
    filenames = [name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]
    number_of_files = len(filenames)
    print(f'Number of isotope-data files in subdirectory: {number_of_files}')
    return number_of_files, filenames, DIR

def file_combiner(nof, filenames, DIR):
    # filenames = filenames[0:4]
    gamma_dict_full, gamma_dict_energy = {},{}
    for fn in range(len(filenames)):
        path = DIR + '//' + filenames[fn]
        df = pd.read_csv(path, header=0)
        df = df[['z', 'n', 'symbol', 'start_level_idx', 'start_level_energy', 'unc_sle', 'end_level_idx', 'end_level_energy', 'unc_ele', 'gamma_idx', 'energy', 'unc_en', 'relative_intensity', 'unc_ri']]
        dfe = df['energy']
        # print(df)
        df_list = df.to_dict(orient='records')
        dfe_list = dfe.to_dict()
        for g in range(len(df)):
            g_key = str(df.at[g,'z']) + '_' + str(df.at[g,'n']) + '_' + str(df.at[g,'start_level_idx']) + '_' + str(df.at[g,'end_level_idx'])
            gamma_dict_full[g_key] = df_list[g]
            gamma_dict_energy[g_key] = dfe_list[g]
        if (g%300 == 0):    
            print(f'Files combined: {round(100*(fn/nof),2)} %')
        
    with open('ENDSF_gamma_energy.json', "w") as json_file:
        json.dump(gamma_dict_energy, json_file, indent=1)
    
    with open('ENDSF_gamma_full.json', "w") as json_file:
        json.dump(gamma_dict_full, json_file, indent=1)
    
    print('Resulting json-files: ENDSF_gamma_energy.json'.ljust(45) + f' | {os.path.getsize('.//ENDSF_gamma_energy.json')/(10**6)} MB'.ljust(12))
    print('Resulting json-files: ENDSF_gamma_full.json'.ljust(45)   + f' | {os.path.getsize('.//ENDSF_gamma_full.json')/(10**6)} MB'.ljust(12))
    

if __name__ == "__main__":
    nof, files, DIR = file_analyser()
    file_combiner(nof, files,DIR)