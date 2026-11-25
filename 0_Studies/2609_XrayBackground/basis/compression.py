import numpy as np


def LL_Compression(data):
    '''
    Dynamic Double Logarithmic Compression (L²C)
    (from Ryan_1988.pdf)
    '''
    
    if type(data) == 'numpy.ndarray':
        pass
    else:
        data = np.array(data)
        
    return np.log(np.log(data+1)+1)

def LLS_Compression(data):
    '''
    Dynamic Double Logarithmic Compression (L²SC)
    (from Morhac_2008.pdf ref 13)
    '''
    
    if type(data) == 'numpy.ndarray':
        pass
    else:
        data = np.array(data)
        
    return np.log(np.log(np.sqrt(data+1)+1)+1)

def LL_Decompression(ddlc_data):
    '''
    Dynamic Double Logarithmic Decompression (L²D)
    Inverse operation to DDLC
    (from Ryan_1988.pdf)
    '''
    
    if type(ddlc_data) == 'numpy.ndarray':
        pass
    else:
        ddlc_data = np.array(ddlc_data)
            
    return np.exp(np.exp(ddlc_data)-1)-1

def LLS_Decompression(ddlc_data):
    '''
    Dynamic Double Logarithmic Compression (L²SD)
    (from Morhac_2008.pdf ref 13)
    '''
    
    if type(ddlc_data) == 'numpy.ndarray':
        pass
    else:
        ddlc_data = np.array(ddlc_data)
        
    return (np.exp(np.exp(ddlc_data)-1)-1)**2 - 1

def ddlc_ddld_error(uncompressed_data, ddlc_data):
    return uncompressed_data-ddlc_data