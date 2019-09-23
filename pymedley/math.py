import numpy as np

################################################################################

def not_normalized_Gaussian(height, half_FWHM, location, x):
        '''

        Purpose:

            Return a Gausian whose height is the variable 'height'

        '''
        FWHM = 2*half_FWHM
        sigma = FWHM / ( 2*np.sqrt( 2*np.log(2) ) )
        gaussian = np.exp( - (x-location)**2 / (2*sigma**2) )
        not_normalized_gaussian = height * gaussian
        return not_normalized_gaussian

def polyfit_with_fixed_points(n, x, y, xf, yf):
    '''

    Purpose:

        Fit a polynomial allowing for fixed points as constraints for the fit

    Taken from:
        
        https://stackoverflow.com/questions/15191088/how-to-do-a-polynomial-fit-with-fixed-points
    
    Author:

      Jaime (Stackoverflow user)
    
    Key word:
      
      Lagrange multipliers

    Useful improvement to do:

        Add one extra Lagrange multiplier to constrain the total area of the fitted polynomial

    '''
    mat = np.empty((n + 1 + len(xf),) * 2)
    vec = np.empty((n + 1 + len(xf),))
    x_n = x**np.arange(2 * n + 1)[:, None]
    yx_n = np.sum(x_n[:n + 1] * y, axis=1)
    x_n = np.sum(x_n, axis=1)
    idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
    mat[:n + 1, :n + 1] = np.take(x_n, idx)
    xf_n = xf**np.arange(n + 1)[:, None]
    mat[:n + 1, n + 1:] = xf_n / 2
    mat[n + 1:, :n + 1] = xf_n.T
    mat[n + 1:, n + 1:] = 0
    vec[:n + 1] = yx_n
    vec[n + 1:] = yf
    params = np.linalg.solve(mat, vec)
    return params[:n + 1][::-1]

def odd_number(number):
    '''Return True is number is odd, False if is even and None if is not integer'''    
    if not isinstance(number, int):
        return None
    else:
        if number%2 == 1:
            return True
        else:
           return False

def even_number(number):
    x = odd_number(number)
    if x is None:
        return None
    else:
        return not x

