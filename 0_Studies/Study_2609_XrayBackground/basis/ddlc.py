import numpy as np


def Dynamic_Double_Logarithmic_Compression(data):
    '''
    Dynamic Double Logarithmic Compression (DDLC)
    (from Ryan_1988.pdf)
    '''
    
    if type(data) == 'numpy.ndarray':
        pass
    else:
        data = np.array(data)
        
    return np.log(np.log(data+1)+1)

def Dynamic_Double_Logarithmic_Decompression(ddlc_data):
    '''
    Dynamic Double Logarithmic Decompression (DDLD)
    Inverse operation to DDLC
    (from Ryan_1988.pdf)
    '''
    
    if type(ddlc_data) == 'numpy.ndarray':
        pass
    else:
        ddlc_data = np.array(ddlc_data)
            
    return np.exp(np.exp(ddlc_data)-1)-1

def ddlc_ddld_error(uncompressed_data, ddlc_data):
    return uncompressed_data-ddlc_data