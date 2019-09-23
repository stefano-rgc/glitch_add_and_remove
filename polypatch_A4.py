#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import tomso.gyre
import tomso.adipls
import tkinter
from IPython import embed
from pymedley.interactive import fit_polynomial

################################################################################

def polypatch_A4(filename, verbose=False):
    '''

    Purpose:
    
        Given an unformated AMDL stellar model, modify the variable A4 within it
        (i.e., the dimensionless N^2) by fitting a polynomial function. The code
        is intended to remove sharp variations in the A4 profile and therefore
        get rid of buoyancy glitches.

        The code has a GUI which allows the user to explore the A4 profile and
        interactively select the interval to fit a polynomial while excluding a
        smaller interval within itself (the latter is suppoused to contain the
        glitch).

        The modification in the A4 profile is done in a consistent way by a
        correspondent sole modification of the first adiabatic exponent, as
        described in the section 2.3 of Ball et al. 2018 ("Surface effects on
        the red giant branch").

    Inputs: 

        - filename (string):

            Name of the unformatted AMDL stellar model file.

    Optional inputs:

        - verbose=False (boolean):

            Print additional information on the terminal. Often useful for
            debugging.
    
    Outputs:

        The program will output a new AMDL stellar model

    Example:

        Let the string 'model.amdl' be the AMDL stellar model

        There are 3 ways to use this program

            1) As an executable on the terminal:

                $ ./polypatch_A4 model.amdl

            2) If you are using a different python kernel from the default one
               (e.g. conda), then it may be better to specify the interpreter
               as, for example:

                $ python3 polypatch_A4 model.amdl

            3) By importing the main function:

                >>>> from polypatch_A4 import polypatch_A4
                >>>> polypatch_A4(model)

    '''

    def window_save():
        '''Ask user whether to save the changes'''

        def save_btn():
            '''Save the modified model into a new AMDL model using the name from the window's entry'''
            filename_output = entry.get()
            tomso.adipls.save_amdl(filename_output, D_modified, A_modified)
            print('Modified AMDL model saved as: {}\n'.format(filename_output))
            window.withdraw()
            window.destroy()

        def close_btn():
            '''Close the tkinter windows without'''
            window.withdraw()
            window.destroy()
            print('Model not saved.\n')

        # Tentative name to save the modified AMDL model
        filename_output_temp = filename+'.modified'

        title = 'Close'
        text = 'Give a name to the modified model. Press cancel to not save.'
        window = tkinter.Tk()
        window.title(title)
        # Text
        tkinter.Label(window, text=text, justify='center').grid(row=0, columnspan=2, sticky='WE', padx=10, pady=10)
        tkinter.Label(window, text='Name').grid(row=1)
        # Entry-text blank field
        entry = tkinter.Entry(window)
        entry.grid(row=1, column=1, sticky='WE')
        entry.insert(0, filename_output_temp)
        # Buttons
        tkinter.Button(window, text='Save', command=save_btn).grid(row=2, column=0, sticky='WE', pady=10, padx=10)
        tkinter.Button(window, text='Cancel', command=close_btn).grid(row=2, column=1, sticky='WE', pady=10, padx=10)
        window.mainloop()

    # Read the AMDL model via the TOMSO package
    D, A = tomso.adipls.load_amdl(filename)
    # Make a copy
    D_modified = D.copy()
    A_modified = A.copy()

    # Ceate the Matplotlib plot for the profile A4 (i.e., the dimensionaless N^2)
    ylabel_adipls = r'$A4\ \equiv\ \frac{r}{g}N^2\ =\ \frac{1}{\Gamma_1}\frac{d \ln p}{d \ln r} - \frac{d \ln \rho}{d \ln r}$'
    xlabel_adipls = r'$A0\ \equiv\ $Fractional radius'
    label_adipls = 'A4'
    fig = plt.figure()
    ax_A4 = plt.axes()
    plot_A4_lines, = ax_A4.plot(A[:,0], A[:,4], color='dodgerblue', linestyle='solid', marker='o', markerfacecolor='dodgerblue', markeredgecolor='None', label=label_adipls)
    ax_A4.grid(True)
    ax_A4.set_yscale('log')
    ax_A4.set_xscale('log')
    ax_A4.set_ylabel(ylabel_adipls, fontsize=20)
    ax_A4.set_xlabel(xlabel_adipls, fontsize=20)
    ax_A4.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)

    # Fit a polynomial in an interval selected interactively
    mask_fit, modified_A = fit_polynomial(fig, ax_A4, plot_A4_lines, verbose=False)

    # Modify the stellar model
    dlnP_dlnr = - ( A[:,3][mask_fit] * A[:,2][mask_fit] ) # Inverse of the pressure heitght scale from Hydrostatic equilibrium
    dlnrho_dlnr = - A[:,2][mask_fit] - A[:,4][mask_fit] # Inverse of the density height scale from the Brun-Vaisala frequency definition
    A_modified[:,4][mask_fit] = modified_A # Mofify the adimensional N2. Make it equal to the fit
    A_modified[:,2][mask_fit] = A[:,2][mask_fit] * A[:,3][mask_fit] # Clean V_g from Gamma_1
    inv_Gamma_1 = (A_modified[:,4][mask_fit] + dlnrho_dlnr) * dlnP_dlnr**-1 # Calculate the reciprocal of the first adiabatic exponent Gamma 1
    A_modified[:,3][mask_fit] = inv_Gamma_1**-1 # Set the new value for Gamma_1 to compensate for the change in N2
    A_modified[:,2][mask_fit] = A_modified[:,2][mask_fit] / A_modified[:,3][mask_fit] # Put back the (now modified) Gamma_1 in V_g

    window_save()

# Called from the terminal (i.e., not imported from another program)
if __name__ == "__main__":
    import sys
    nargs = len(sys.argv)
    if nargs == 2:
        stellar_model = sys.argv[1]
        polypatch_A4(stellar_model, verbose=False)
    else:
        print('! {} arguments given. Please give the name of the AMDL stellar model as only argument.\n'.format(nargs-1))