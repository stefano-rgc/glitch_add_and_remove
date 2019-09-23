import pandas as pd  
import numpy as np    

def read_famdl(filename):
    '''
    Then the radius is famdl[:,0]
    the Gamma1 is famdl[:,3]
    etc
    (as in CoRoT_ESTA_Files.pdf section FAMDL)
    '''
                                                
    famdl = pd.read_fwf(filename, skiprows=2, widths=[20,20,20,20])    
    famdl = famdl.values.reshape(-1)                                                                                                              
    famdl = famdl[~np.isnan(famdl)].reshape(-1,6)
    return famdl
