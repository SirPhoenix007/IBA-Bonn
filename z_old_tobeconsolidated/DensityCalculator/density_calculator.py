import numpy as np
from mendeleev import element


def density_calc(composition: dict):
    elem = list(composition.keys())
    percentages = list(composition.values())
    result_density = 0

    print(elem)
    
    for i in range(len(elem)):
        result_density += element(elem[i]).density * (percentages[i]/100)
    
    print(result_density)
    
    return result_density

blade15front = {"Cu":92.6, "As":2.20, "Ag":0.43, "Sn":0.68, "Sb":0.43, "Pb":0.58, "Ni":2.90, "Zn":0.15}
blade15back = {"Cu":94.2, "As":1.98, "Ag":0.37, "Sn":0.61, "Sb":0.19, "Pb":0.26, "Ni":2.40}

density_calc(blade15front)

density_calc(blade15back)