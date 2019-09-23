#!/usr/bin/env python3
import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets 
import tomso.gyre
import tomso.adipls
from IPython import embed
import tkinter
from tkinter import Tk, simpledialog
from tkinter.filedialog import askopenfilename
import pickle
from copy import deepcopy
from pymedley.interactive import ask_output_name
from pymedley.math import odd_number
from pymedley.mpl import clear_line2D, axvlines
from pymedley.meta import get_kwargs
from pymedley.interactive import Textbox, Button, Checkbox
from pymedley.astero import Read_stellar_model

class Make_plots:
    '''Object containing the plots from the instance of the class Read_stellar_model'''

    def __init__(self, stellar_model):
        if not stellar_model.type == 'adipls':
            raise ValueError("Support for ADIPLS model only. Your model is '{}'.".format(stellar_model.type))
        self.stellar_model_type = stellar_model.type
        # Create two subplots: N^2 VS the normalized buoyancy depth  and for A[:,4] VS fractional radius
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 6))
        self.figure = fig
        self.ax1 = ax1
        self.ax2 = ax2
        self.make_plots(stellar_model, Verbose=False, from_init=True)

    def make_plots(self,\
                   stellar_model,\
                   xlabel_fontsize=15,\
                   ylabel_fontsize=20,\
                   legend_fontsize=10,\
                   grid=True,\
                   scale='log',\
                   ylabel_ax1=r'$N^2$ (Hz$^2$)',\
                   xlabel_ax1=r'Buoyancy depth $(\omega_g)$',\
                   label_ax1=r'$N^2$',\
                   label_ax2_adipls=r'$A4$',\
                   ylabel_ax2_adipls=r'$A4\equiv\frac{r}{g}N^2$',\
                   Verbose=False,
                   from_init=False):
        '''Make plots'''

        if stellar_model.type != self.stellar_model_type:
            raise ValueError('The type of stellar model from the instance Read_stellar_model does not match the one from the given argument.')

        try: # Clear existing curves
            clear_line2D(self.figure, self.ax1_line2D, self.ax1)
            clear_line2D(self.figure, self.ax2_line2D, self.ax2)
        except AttributeError as e:
            if Verbose:
                print('Harmless AttributeError :::', e)

        # Left plot: normalized buoyancy depth VS N^2. Limited to the buoyancy cavity
        plot_line2D1, = self.ax1.plot(stellar_model.normalized_buoyancy_depth,\
                                      stellar_model.N2_buoyancy,\
                                      color='dodgerblue', linestyle='solid', marker='o',\
                                      markerfacecolor='dodgerblue', markeredgecolor='None',\
                                      label=label_ax1)

        # Right plot: fractional_r VS A4. Limited to the buoyancy cavity
        if self.stellar_model_type == 'adipls':
            self.ax2.set_ylabel(ylabel_ax2_adipls, fontsize=ylabel_fontsize)
            plot_line2D2, = self.ax2.plot(stellar_model.fractional_r_buoyancy,\
                                          stellar_model.A4_buoyancy,\
                                          color='dodgerblue', linestyle='solid', marker='o',\
                                          markerfacecolor='dodgerblue', markeredgecolor='None',\
                                          label=label_ax2_adipls)

        self.ax1.set_ylabel(ylabel_ax1, fontsize=ylabel_fontsize)
        self.ax1.set_xlabel(xlabel_ax1, fontsize=xlabel_fontsize)
        self.ax1.legend(loc='best', ncol=1, framealpha=0.5, fontsize=legend_fontsize)

        if from_init:
            self.ax1.invert_xaxis()
            self.ax1.grid(grid)
            self.ax1.set_yscale(scale)
            self.ax1.set_xscale(scale)
            self.ax2.grid(grid)
            self.ax2.set_yscale(scale)
            self.ax2.set_xscale(scale)
        else:
            self.ax1.set_yscale(self.ax1.get_yscale())
            self.ax1.set_xscale(self.ax1.get_xscale())
            self.ax2.set_yscale(self.ax2.get_yscale())
            self.ax2.set_xscale(self.ax2.get_xscale())

        self.ax2.set_xlabel(r'$A0 \equiv r/R$',fontsize=xlabel_fontsize)
        self.ax2.legend(loc='best', ncol=1, framealpha=0.5, fontsize=legend_fontsize)

        # Collect specific information
        self.ax1_line2D = plot_line2D1
        self.ax2_line2D = plot_line2D2

    def redraw(self):
        self.ax1.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)
        self.ax2.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)
        self.figure.canvas.draw()

class Gaussian_parameters:
    ''' '''
    def __init__(self, loc, sigma, hat_A_G, name=''):
        ''' '''
        self.name = name

        # Direct parameters
        self.loc = loc
        self.sigma = sigma
        self.hat_A_G = hat_A_G

        # Derivated parameters
        self.amplitude = self.hat_A_G / ( np.sqrt(2*np.pi) * self.sigma )
        self.ratio_N2_N02 = self.amplitude + 1

    def update_loc(self, loc):
        self.loc = loc

    def update_sigma(self, sigma):
        self.sigma = sigma
        # Consequence
        self.amplitude = self.hat_A_G / ( np.sqrt(2*np.pi) * sigma )

    def update_hat_A_G(self, hat_A_G):
        self.hat_A_G = hat_A_G
        # Consequence
        self.amplitude = hat_A_G / ( np.sqrt(2*np.pi) * self.sigma )

    def print(self,\
              tabs=0,\
              comment='',\
              text_hat_A_G='hat A_G',\
              text_sigma='Sigma',\
              text_loc='Location',\
              hat_A_G=True,\
              sigma=True,\
              loc=True):

        tabs =  '\t' * tabs
        print( tabs + 'Gaussian parameters:' )
        if self.name != '': 
            print( tabs + 'Name: {}'.format(self.name) )
        if comment != '':
            print ( tabs    + 'Comment: "{}"'.format(comment) )
        if loc:
            print( tabs + '{} = {}'.format(text_loc, self.loc) )        
        if sigma:
            print( tabs + '{} = {}'.format(text_sigma, self.sigma) )        
        if hat_A_G:
            print( tabs + '{} = {}'.format(text_hat_A_G, self.hat_A_G) )

class Interactive_variables:
    '''Organize the interactive variables'''
    keystroke = None
    pressed_button = None
    released_button = None
    pressed_xdata = None
    pressed_ydata = None
    released_xdata = None
    released_ydata = None
    wrt_final_buoyancy_depth = False
    new_resolution = False

def update_kwarg(funct, kwargs_new):
    '''
    Given a dictionary and a function, if any of the keys of the former matches
    the name of a kwarg of the latter, give the value of that key to that kwarg.
    '''
    kwargs = get_kwargs(funct)
    for key in kwargs:
        if key in kwargs_new:
            kwargs[key] = kwargs_new[key]
    return kwargs

def scaled_Gaussian(x, height, loc, sigma):
    '''Return a scaled Gausian'''
    gaussian = np.exp( - (x-loc)**2 / (2*sigma**2) ) # Gaussian with out the normalization factor (maximum value is 1)
    scaled_gaussian = height * gaussian
    return scaled_gaussian

def parse_user_command_line_input(args):
    '''Parse the command line. Identify arguments'''
    
    # Meaning of the one-letter flags
    available_flags = { '-t':'stellar_model_type',\
                        '-b':'buoyancy_cavity_between',\
                        '-s':'save',\
                        '-v':'Verbose',\
                        '-G':'G',\
                        '-p':'gaussian_peak_initial_text',\
                        '-A':'gaussian_hat_A_G_initial_text',\
                        '-w':'gaussian_sigma_initial_text',\
                        '-l':'gaussian_loc_initial_text',\
                        '-i':'interactive',\
                        '-o':'filename_output',\
                        '-f':'wrt_final_buoyancy_depth',\
                        '-r':'extra_meshpoints'}

    # Pre-defined posible values input by the user
    available_flag_values = { '-t':{'adipls', 'amdl', 'gyre'},\
                              '-s':{'n','no','not','y','yes','true','false'},\
                              '-v':{'n','no','not','y','yes','true','false'},\
                              '-i':{'n','no','not','y','yes','true','false'},\
                              '-f':{'n','no','not','y','yes','true','false'} }

    # Conversion of the values input by the user to the actual values
    available_flag_values_dictionary = { 'yes':True,\
                                         'y':True,\
                                         'true':True,\
                                         'not':False,\
                                         'no':False,\
                                         'n':False,\
                                         'false':False }
    
    # Check for the correct number of arguments
    narg = len(args)
    if odd_number(narg):
        raise ValueError("**Wrong number of command line inputs. Please see the help.**\n")

    # Check if help is requested
    if narg == 2 and args[1] == '-h':
        print(add_gaussian_A4.__doc__)
        raise SystemExit(0) 

    # Get the filename
    filename = str(args[1])

    # Assume the flags start form the third comand-line argument
    flags = args[2::2]
    flags_value = args[3::2]
    flags_value = [ str(fv).casefold() for fv in flags_value ] # Make lowercase
    
    # Check for unknown flags
    if not set(flags).issubset( available_flags.keys() ):
        raise ValueError("**Wrong flag as input. Please see the help.**")
    
    # Check for repeated flags
    if len( set(flags) ) < len(flags):
        raise ValueError("**Repeated flags.**\n")

    # Check for valid flag values
    for f,fv in zip(flags, flags_value):
        # Exclude flags without default values like -G or -l
        if f in available_flag_values:
            if not fv in available_flag_values[f]:
                raise ValueError( "**Please, select a valid input among {} for the flag '{}'.**".format( available_flag_values[f], f ) )

    # Collect flags and values
    flags_and_values = { f:v for f,v in zip(flags,flags_value) }
    
    # Check if the type of stellar model is set            
    if '-t' in flags:
        stellar_model_type = flags_and_values.pop('-t')
    else:
        # Figure out the type of tellar model
        if re.match(r'.*gyre.*',filename, re.IGNORECASE):
            stellar_model_type = 'gyre'
            print( "**The stellar model '{}' is assumed to have the MESA-GYRE format.**".format(filename) )
        elif re.match(r'.*amdl.*',filename, re.IGNORECASE):
            stellar_model_type = 'amdl'
            print( "**The stellar model '{}' is assumed to have the MESA-GYRE format.**".format(filename) )
        elif re.match(r'.*adipls.*',filename, re.IGNORECASE):
            stellar_model_type = 'adipls'
            print( "**The stellar model '{}' is assumed to be the amdl unfformated file.**".format(filename) )
        else:
            raise ValueError("**Please, specify the type of stellar model (amdl or gyre) by the flag '-t'.**")

    # Limits of the interval to search for the buoyancy cavity (units of fractional radius)
    if not '-b' in flags:
        raise ValueError("**The flag '-b' is mandatory. Please, specify the interval to search for the buoyancy cavity in units of fractional radius.**")
    # Extract the flag from the dictionary
    temp = flags_and_values.pop('-b')
    temp = temp.split(',')
    # Check it is two numbers separated by a comma
    if not len(temp)==2:
        raise ValueError("**Input two valid numbers in units of fractional radius separated by a comma and no space to bound the buoyancy cavity.**")
    # Convert the string into number
    try:
        buoyancy_cavity_between = ( np.float(temp[0]), np.float(temp[1]) )
    except ValueError:
        raise ValueError("**Input two valid numbers in units of fractional radius separated by a comma and no space to search for the buoyancy cavity.**")
        
    # Check if a value for increment the number of meshpoints is set
    if '-r' in flags:
        try:
            flags_and_values['-r'] = np.int( flags_and_values['-r'] )
        except ValueError:
            raise ValueError("**Input a valid integer for the number of extra meshpoints in between the original mesh. This will apply for a range of 3 sigma about the center of the added Gaussian.**")

    # Check if a value for the gravitational constant is set
    if '-G' in flags:
        try:
            flags_and_values['-G'] = np.float( flags_and_values['-G'] )
        except ValueError:
            raise ValueError("**Input a valid number for the gravitational constant in cgs units (e.g. 6.67428e-08)**")

    # Check if the location for the Gaussian-like function is set
    if '-l' in flags:
        try:
            flags_and_values['-l'] = np.float( flags_and_values['-l'] )
        except ValueError:
            raise ValueError("**Input a valid number for the location of the Gaussian-like function in units of normalized buoyancy depth**")

    # Check if the value of the peak of the Gaussian-like function is set
    if '-p' in flags:
        try:
            flags_and_values['-p'] = np.float( flags_and_values['-p'] )
        except ValueError:
            raise ValueError("**Input a valid number for the value of the peak of the Gaussian-like function in units of Hz^2**")

    # Check if the value of the half of the FWHM of the Gaussian-like function is set
    if '-w' in flags:
        try:
            flags_and_values['-w'] = np.float( flags_and_values['-w'] )
        except ValueError:
            raise ValueError("**Input a valid number for the half of the FWHM of the Gaussian-like function in units of normalized buoyancy depth**")

    # Check if the parameter for the amplitud of the Gaussian-like function is set (as described in Cunha et. al.'s 2019 paper)
    if '-A' in flags:
        try:
            flags_and_values['-A'] = np.float( flags_and_values['-A'] )
        except ValueError:
            raise ValueError("**Input a valid number for the amplitud of the Gaussian-like function in units of normalized buoyancy depth (as described in Cunha et. al.'s 2019 paper)**")

    # Convert the user command-line flag-value input into proper Python-language variables
    for f in ['-s', '-v', '-i', '-f']:
        if f in flags_and_values:
            flags_and_values[f] = available_flag_values_dictionary[ flags_and_values[f] ]

    # If the interactive mode is off, check that loc, sigma and amplitude of the gaussian are set
    if '-i' in flags_and_values:
        if flags_and_values['-i']==False:
            if not '-l' in flags_and_values:
                raise ValueError("**Because the interactive mode is off, please specify a location for the Gaussian using the flag '-l'**")
            if not '-w' in flags_and_values:
                raise ValueError("**Because the interactive mode is off, please specify a sigma for the Gaussian using the flag '-w'**")
            if not '-A' in flags_and_values:
                raise ValueError("**Because the interactive mode is off, please specify an amplitude for the Gaussian using the flag '-A'**")

    # Convert the user command-line flag input accordingly to available_flags (e.g., '-l' to 'gaussian_loc_initial_text')
    flags_and_values = { available_flags[f]:fv for f,fv in flags_and_values.items() }

    # Return parsed command line
    return filename, stellar_model_type, buoyancy_cavity_between, flags_and_values                    

def glitch_wrt_model(gaussian_loc, stellar_model, Verbose=False, comment=''):
    '''Compute values of the stellar model at the position od the glitch'''
    
    if np.all(np.diff(stellar_model.normalized_buoyancy_depth) > 0):
        r_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth, stellar_model.r_buoyancy)
        N_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth, stellar_model.N_buoyancy)
        N2_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth, stellar_model.N2_buoyancy)
    else:
        r_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth[::-1], stellar_model.r_buoyancy[::-1])
        N_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth[::-1], stellar_model.N_buoyancy[::-1])
        N2_star = np.interp(gaussian_loc, stellar_model.normalized_buoyancy_depth[::-1], stellar_model.N2_buoyancy[::-1])

    modulation_r_buoyancy = np.sqrt( stellar_model.r_buoyancy / r_star )
    modulation_N_buoyancy = np.sqrt( N_star / stellar_model.N_buoyancy )
    gaussian_modulation_buoyancy =  modulation_r_buoyancy * modulation_N_buoyancy

    # Calculate the value of the modulation at the glitch position in order to plot the modulation normalized to one at the glitch position
    modulation_star = np.interp(r_star, stellar_model.r_buoyancy, gaussian_modulation_buoyancy)

    if Verbose:
        # Print starred valuer
        print('')
        print('[from glitch_wrt_model():]')
        print('|Starred values:' )
        if comment != '':
            print('|Comment: ' + comment )
        print('|r = {:.3}'.format(r_star) )
        print('|N = {:.3}'.format(N_star) )
        print('|N2 = {:.3}'.format(N2_star) )
        print('')

    return { 'r':r_star,\
             'N':N_star,\
             'N2':N2_star,\
             'modulation_r_buoyancy':modulation_r_buoyancy,\
             'modulation_N_buoyancy':modulation_N_buoyancy,\
             'gaussian_modulation_buoyancy':gaussian_modulation_buoyancy,\
             'modulation_starred':modulation_star }

def add_gaussian_interactively(stellar_model,\
                               plots,\
                               gaussian_peak_initial_text='',\
                               gaussian_hat_A_G_initial_text=0.0,\
                               gaussian_sigma_initial_text=0.0,\
                               gaussian_loc_initial_text=0.0,\
                               Verbose=True,\
                               interactive=True,\
                               wrt_final_buoyancy_depth=False,\
                               extra_meshpoints=None):
    '''

    Purpose:

        Define the parameters of the Gaussian-like buoyancy-glitch to add to the
        stellar model

    Inputs:

        stellar_model

            An instance if the class Read_Stellar_model

        plots

            An instance of the class Make_plots

    '''

    def onPressKey(event):
        ''' Get the key pressed during the plot visualization '''
        iinput.keystroke = event.key
        if Verbose:
            print('\nUpdate of variable (by pressing key): "iinput.keystroke" <--- ', iinput.keystroke, '\n')

    def onPressButton(event):
        ''' Get which mouse button is pressed and the xdata of the location where the mouse button is pressed during the plot visualization '''
        iinput.pressed_button = event.button
        iinput.pressed_xdata = event.xdata
        iinput.pressed_ydata = event.ydata
        if Verbose:
            print('')
            print('Update of variable (by mouse): "iinput.pressed_button" <--- ', iinput.pressed_button)
            print('Update of variable (by mouse): "iinput.pressed_xdata" <--- ', iinput.pressed_xdata)
            print('Update of variable (by mouse): "iinput.pressed_ydata" <--- ', iinput.pressed_ydata)
            print('')

    def onReleaseButton(event):
        ''' Get the ydata of the location where the mouse button is released during the plot visualization '''
        iinput.released_button = event.button
        iinput.released_xdata = event.xdata
        iinput.released_ydata = event.ydata
        if Verbose:
            print('')
            print('Update of variable (by mouse): "iinput.released_button" <--- ', iinput.released_button)
            print('Update of variable (by mouse): "iinput.released_xdata" <--- ', iinput.released_xdata)
            print('Update of variable (by mouse): "iinput.released_ydata" <--- ', iinput.released_ydata)
            print('')    

    def read_peak(text):
        '''Get the value for the peak of N^2 from the text box'''
        itextbox['N2peak'] = text
        if text == '':
            results['gaussian_N2'] = np.nan
        else:
            try:
                results['gaussian_N2'] = np.float(text)
                textbox_peak.textbox.color='0.95'
            except ValueError:
                print( '\t' + '*Input a valid number for the Gaussian peak*' )
                results['gaussian_N2'] = np.nan
                textbox_peak.textbox.color='red'
        if Verbose:
            print('itextbox[N2peak] <---', text)
            print("results['gaussian_N2'] = ",  results['gaussian_N2'])

    def read_sigma(text):
        '''Get the value for the sigma from the text box'''
        itextbox['sigma'] = text
        if text == '':
            gaussian_p.update_sigma(np.nan)
        else:
            try:
                gaussian_p.update_sigma( np.float(text) )
                textbox_sigma.textbox.color='0.95'
            except ValueError:
                print( '\t' + '*Input a valid number for the Gaussian sigma*' )
                gaussian_p.update_sigma(np.nan)
                textbox_sigma.textbox.color='red'
        if Verbose:
            print('itextbox[sigma] <---', text)
            print("gaussian_p.sigma = ",  gaussian_p.sigma)

    def read_loc(text):
        '''Get the value for the location from the text box'''
        itextbox['loc'] = text
        if text == '':
            gaussian_p.update_loc(np.nan)
        else:
            try:
                gaussian_p.update_loc( np.float(text) )
                textbox_loc.textbox.color='0.95'
            except ValueError:
                print( '\t' + '*Input a valid number for the Gaussian location*' )
                gaussian_p.update_loc(np.nan)
                textbox_loc.textbox.color='red' 
        if Verbose:
            print('itextbox[loc] <---', text)
            print("gaussian_p.loc = ",  gaussian_p.loc)

    def read_hat_A_G(text):
        '''Get the value for the amplitude from the text box'''
        itextbox['hat_A_G'] = text
        if text == '':
            gaussian_p.update_hat_A_G(np.nan)
        else:
            try:
                gaussian_p.update_hat_A_G( np.float(text) )
                textbox_hat_A_G.textbox.color='0.95'
            except ValueError:
                print( '\t' + '*Input a valid number for the Gaussian hat A_G*' )
                gaussian_p.update_hat_A_G(np.nan)
                textbox_hat_A_G.textbox.color='red'
        if Verbose:
            print('itextbox[hat_A_G] <---', text)
            print("gaussian_p.hat_A_G = ",  gaussian_p.hat_A_G)

    def wrt_final_buoyancy_depth_check(label):
        '''Check if the user has marked the option w.r.t. the final buoyancy depth'''
        iinput.wrt_final_buoyancy_depth = checkbox_wrt_final.checkbox.get_status()[0]
        if Verbose:
            print("iinput.wrt_final_buoyancy_depth = ", iinput.wrt_final_buoyancy_depth)
        plots.redraw()

    def add_gaussian(event):
        '''Add the Gaussian-like glitch to the N^2 profile'''

        # Check flags
        wrt_final_buoyancy_depth_flag = wrt_final_buoyancy_depth
        new_resolution_flag = iinput.new_resolution if interactive else bool(extra_meshpoints)

        # Check that the interval where to draw the Gaussian is selected
        if not gaussian_interval['ax1']['ready']:
            print('*Please select an interval before adding a Gaussian*')
        else:
            # For convinience, give an alias to the mask containing the Gaussian
            mask_gaussian = gaussian_interval['ax1']['mask']

            # CASE (1): The peak in N^2 of the new Gaussian is not given
            if itextbox['N2peak'] is ''\
               and np.isfinite(gaussian_p.loc)\
               and np.isfinite(gaussian_p.sigma)\
               and np.isfinite(gaussian_p.hat_A_G):
                starred_values = glitch_wrt_model(gaussian_p.loc, stellar_model.new_resolution if new_resolution_flag else stellar_model)
                record_case = 'Output parameters for the Gaussian [inputs parameters: loc, sigma and hat A_G]:'

            # CASE (2): The sigma of the new Gaussian is not given
            elif np.isnan(gaussian_p.sigma)\
                 and np.isfinite(gaussian_p.loc)\
                 and np.isfinite(gaussian_p.hat_A_G)\
                 and np.isfinite(np.float(itextbox['N2peak'])):
                starred_values = glitch_wrt_model(gaussian_p.loc, stellar_model.new_resolution if new_resolution_flag else stellar_model)
                # Calculate the value of the peak in N2 (at the position gaussian_p.loc)
                results['gaussian_N2'] = np.float( itextbox['N2peak'] )
                # Ratio between N2 value with and without gaussian at the position gaussian_parameter.loc 
                ratio_N2_N02 =  results['gaussian_N2'] / starred_values['N2']
                # Calculate sigma
                gaussian_sigma = ( 1 / np.sqrt(2*np.pi) ) * ( gaussian_p.hat_A_G / (ratio_N2_N02 - 1) )
                gaussian_p.update_sigma(gaussian_sigma)
                record_case = 'Output parameters for the Gaussian [inputs papameters: loc, peank N2 and hat A_G]:'

            # CASE (3): The amplitude of the new Gaussian is not given
            elif np.isnan(gaussian_p.hat_A_G)\
                 and np.isfinite(gaussian_p.loc)\
                 and np.isfinite(gaussian_p.sigma)\
                 and np.isfinite(np.float(itextbox['N2peak'])):
                starred_values = glitch_wrt_model(gaussian_p.loc, stellar_model.new_resolution if new_resolution_flag else stellar_model)
                # Calculate the value of the peak in N2 (at the position gaussian_p.loc)
                results['gaussian_N2'] = np.float( itextbox['N2peak'] )
                # Ratio between N2 value with and without gaussian at the position gaussian_parameter.loc 
                ratio_N2_N02 =  results['gaussian_N2'] / starred_values['N2']
                # Calculate hat A_G
                gaussian_hat_A_G = np.sqrt(2*np.pi) * gaussian_p.sigma * (ratio_N2_N02 - 1)
                gaussian_p.update_hat_A_G(gaussian_hat_A_G)
                record_case = 'Output parameters for the Gaussian [input parameters: loc, sigma and peak N2]:'

            # CASE (4): 
            else:
                raise ValueError('**The peak in N^2, sigma and hat A_G of the Gaussian parameters cannot be all set at the same time. Please correct it.**')

            # Once the Gaussian parameters are set, generate the Gaussian
            results['new_N2_gaussian'] = get_new_N2(mask_gaussian=mask_gaussian,\
                                                    stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model,\
                                                    gaussian_modulation_buoyancy=starred_values['gaussian_modulation_buoyancy'],\
                                                    gaussian_parameters=gaussian_p,\
                                                    comment=record_case)

            if interactive:
                wrt_final_buoyancy_depth_flag = iinput.wrt_final_buoyancy_depth

            if wrt_final_buoyancy_depth_flag:
                
                # Initialize variables for the do-while loop
                iteration = 1
                max_iteration = 20
                gaussian_p_perturbed_model_aimed = deepcopy(gaussian_p)

                # Convergence level
                convergence_value_gaussian_loc     = 0.0000001
                convergence_value_gaussian_sigma   = 0.0000001
                convergence_value_gaussian_hat_A_G = 0.0000001

                # Steps to compute the derivatives w.r.t. the Gaussian loc, sigma and hat A_G
                h_gaussian_loc     = 0.00001 
                h_gaussian_sigma   = 0.00001
                h_gaussian_hat_A_G = 0.00001

                # do-while loop
                while True:

                    # Get the Gaussian parameters w.r.t. the perturbed model
                    gaussian_parameters_perturbed_model = gaussian_parameters_wrt_final_model(stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model,\
                                                                                              new_N2_gaussian=results['new_N2_gaussian'],\
                                                                                              gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                              starred_values_wrt_original=starred_values,\
                                                                                              plot=False)

                    # Compare this value to the aimed one
                    difference_gaussian_loc     =  gaussian_parameters_perturbed_model['loc']     - gaussian_p_perturbed_model_aimed.loc
                    difference_gaussian_sigma   =  gaussian_parameters_perturbed_model['sigma']   - gaussian_p_perturbed_model_aimed.sigma
                    difference_gaussian_hat_A_G =  gaussian_parameters_perturbed_model['hat_A_G'] - gaussian_p_perturbed_model_aimed.hat_A_G
                    relative_difference_gaussian_loc     = difference_gaussian_loc     / gaussian_p_perturbed_model_aimed.loc   
                    relative_difference_gaussian_sigma   = difference_gaussian_sigma   / gaussian_p_perturbed_model_aimed.sigma   
                    relative_difference_gaussian_hat_A_G = difference_gaussian_hat_A_G / gaussian_p_perturbed_model_aimed.hat_A_G   
                    print('Iteration:', iteration)                       
                    print('(!) Relative difference Gaussian location = ', relative_difference_gaussian_loc)
                    print('(!) Relative difference Gaussian sigma = ', relative_difference_gaussian_sigma)
                    print('(!) Relative difference Gaussian hat A_G = ', relative_difference_gaussian_hat_A_G)

                    # Assess convergene
                    c1 = np.abs(relative_difference_gaussian_loc)     < convergence_value_gaussian_loc
                    c2 = np.abs(relative_difference_gaussian_sigma)   < convergence_value_gaussian_sigma
                    c3 = np.abs(relative_difference_gaussian_hat_A_G) < convergence_value_gaussian_hat_A_G
                    c4 = iteration > max_iteration
                    break_condition = (c1 and c2 and c3) or c4
                    if break_condition:
                        break

                    # If while continues, then prepare requirements for computing next Gaussian-parameter values with the Newton-Raphson method

                    # Obtain the derivatives: How much the perturbed-model  Gaussian-parameters change when changing the nperturbed-model Gaussian-parameters
                    derivatives = get_derivatives(unperturbed_stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model,\
                                                  gaussian_parameters_0=gaussian_p,\
                                                  h_loc=h_gaussian_loc,\
                                                  h_sigma=h_gaussian_sigma,\
                                                  h_hat_A_G=h_gaussian_hat_A_G,\
                                                  mask_gaussian=mask_gaussian)

                    # Jacovian
                    jacobian = np.array([ [ derivatives['wrt_loc']['loc'],     derivatives['wrt_sigma']['loc'],     derivatives['wrt_hat_A_G']['loc']     ] ,\
                                          [ derivatives['wrt_loc']['sigma'],   derivatives['wrt_sigma']['sigma'],   derivatives['wrt_hat_A_G']['sigma']   ] ,\
                                          [ derivatives['wrt_loc']['hat_A_G'], derivatives['wrt_sigma']['hat_A_G'], derivatives['wrt_hat_A_G']['hat_A_G'] ] ])

                    # Inverse of Jacovian
                    jacobian_inv = np.linalg.inv(jacobian)
                    
                    # Pivor value vector
                    gaussian_parameters_vector = np.array([ [gaussian_p.loc]     ,\
                                                            [gaussian_p.sigma]   ,\
                                                            [gaussian_p.hat_A_G] ])
                    
                    # Vectorial function to be found the root
                    differences_vector = np.array([ [difference_gaussian_loc]     ,\
                                                    [difference_gaussian_sigma]   ,\
                                                    [difference_gaussian_hat_A_G] ])

                    # Newton Raphson method
                    new_gaussian_parameters_vector = gaussian_parameters_vector - np.matmul( jacobian_inv, differences_vector )

                    # Update the gaussian parameters by increassin loc by h
                    gaussian_p.update_loc(new_gaussian_parameters_vector[0][0])
                    gaussian_p.update_sigma(new_gaussian_parameters_vector[1][0])
                    gaussian_p.update_hat_A_G(new_gaussian_parameters_vector[2][0])

                    starred_values = glitch_wrt_model(gaussian_p.loc,stellar_model.new_resolution if new_resolution_flag else stellar_model)

                    # Compute again the new N2
                    results['new_N2_gaussian'] = get_new_N2(mask_gaussian=mask_gaussian,\
                                                            stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model,\
                                                            gaussian_modulation_buoyancy=starred_values['gaussian_modulation_buoyancy'],\
                                                            gaussian_parameters=gaussian_p,\
                                                            comment=record_case)

                    # Record iteration
                    iteration += 1


            ax2x, ax2y = get_data_for_ax2(mask_gaussian=mask_gaussian,\
                                          stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model)

            # Calculate the value of the peak in N2 (at the position gaussian_p.loc)
            # !!!!!!!*** Make sure this does not contradict the case when the peak N2 is given ***!!!!!!!
            if np.all(np.diff(stellar_model.normalized_buoyancy_depth) > 0):
                results['gaussian_N2'] = np.interp( gaussian_p.loc, stellar_model.new_resolution.normalized_buoyancy_depth[mask_gaussian] if new_resolution_flag else stellar_model.normalized_buoyancy_depth[mask_gaussian], results['new_N2_gaussian'] )
            else:
                results['gaussian_N2'] = np.interp( gaussian_p.loc, stellar_model.new_resolution.normalized_buoyancy_depth[mask_gaussian][::-1] if new_resolution_flag else stellar_model.normalized_buoyancy_depth[mask_gaussian][::-1], results['new_N2_gaussian'][::-1] ) 

            # If interactive version, then do and show the plots
            if interactive:
                plot_gaussian(mask_gaussian=mask_gaussian,\
                              normalized_buoyancy_depth=stellar_model.new_resolution.normalized_buoyancy_depth if new_resolution_flag else stellar_model.normalized_buoyancy_depth,\
                              gaussian_modulation_buoyancy=starred_values['gaussian_modulation_buoyancy'],\
                              modulation_scale=results['gaussian_N2']/starred_values['modulation_starred'],\
                              ax2x=ax2x,\
                              ax2y=ax2y)

            results['ax2y'] = ax2y

            gaussian_parameters_perturbed_model= gaussian_parameters_wrt_final_model(stellar_model=stellar_model.new_resolution if new_resolution_flag else stellar_model,\
                                                                                     new_N2_gaussian=results['new_N2_gaussian'],\
                                                                                     gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                     starred_values_wrt_original=starred_values,\
                                                                                     plot=False) # Set plot=True to test if wrt_final_model is sensible

    def get_derivatives(unperturbed_stellar_model,\
                        gaussian_parameters_0,\
                        h_loc,\
                        h_sigma,\
                        h_hat_A_G,\
                        mask_gaussian):
        ''' Compute how much the quantities in the perturbed model chanche due to a different position of the gaussian in the unperturbated model '''

        # Gaussian parameters for which the derivate is desired
        #gaussian_parameters_0 = Gaussian_parameters(loc=loc, sigma=sigma, hat_A_G=hat_A_G)
        # Get the starred values (i.e., stellar model's values at the position of the glitch)
        starred_values_0 = glitch_wrt_model(gaussian_parameters_0.loc,unperturbed_stellar_model)
        # Create the Gaussian for the Gaussian parameters *without* increment
        N2_gaussian_0 = get_new_N2(mask_gaussian=mask_gaussian,\
                                   stellar_model=unperturbed_stellar_model,\
                                   gaussian_modulation_buoyancy=starred_values_0['gaussian_modulation_buoyancy'],\
                                   gaussian_parameters=gaussian_parameters_0)
        # Obtain the Gaussian parameters w.r.t. the perturbed model
        gaussian_parameters_perturbed_model_0 = gaussian_parameters_wrt_final_model(stellar_model=unperturbed_stellar_model,\
                                                                                    new_N2_gaussian=N2_gaussian_0,\
                                                                                    gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                    starred_values_wrt_original=starred_values_0,\
                                                                                    plot=False)

        ### W.r.t gaussian loc ### ---------------------------------------------------------

        gaussian_parameters_h_loc = Gaussian_parameters(loc=gaussian_parameters_0.loc+h_loc, sigma=gaussian_parameters_0.sigma, hat_A_G=gaussian_parameters_0.hat_A_G)
        # Get the starred values (i.e., stellar model's values at the position of the glitch)
        starred_values_h_loc = glitch_wrt_model(gaussian_parameters_h_loc.loc,unperturbed_stellar_model)
        # Create the Gaussian for the Gaussian parameters *with* increment
        N2_gaussian_h_loc = get_new_N2(mask_gaussian=mask_gaussian,\
                                       stellar_model=unperturbed_stellar_model,\
                                       gaussian_modulation_buoyancy=starred_values_h_loc['gaussian_modulation_buoyancy'],\
                                       gaussian_parameters=gaussian_parameters_h_loc)
        # Obtain the Gaussian parameters w.r.t. the perturbed model
        gaussian_parameters_perturbed_model_h_loc = gaussian_parameters_wrt_final_model(stellar_model=unperturbed_stellar_model,\
                                                                                        new_N2_gaussian=N2_gaussian_h_loc,\
                                                                                        gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                        starred_values_wrt_original=starred_values_h_loc,\
                                                                                        plot=False)

        # Compute the derivative of the Gaussian parameters w.r.t. the Gaussian's location 
        dhat_A_G_dloc = ( gaussian_parameters_perturbed_model_h_loc['hat_A_G'] - gaussian_parameters_perturbed_model_0['hat_A_G'] ) / h_loc
        dsigma_dloc   = ( gaussian_parameters_perturbed_model_h_loc['sigma']   - gaussian_parameters_perturbed_model_0['sigma'] )   / h_loc
        dloc_dloc     = ( gaussian_parameters_perturbed_model_h_loc['loc']     - gaussian_parameters_perturbed_model_0['loc']   )   / h_loc

        # Collect the derivatives
        derivatives_wrt_loc = { 'hat_A_G':dhat_A_G_dloc, 'sigma':dsigma_dloc, 'loc':dloc_dloc }

        ### W.r.t gaussian sigma ### ---------------------------------------------------------

        gaussian_parameters_h_sigma = Gaussian_parameters(loc=gaussian_parameters_0.loc, sigma=gaussian_parameters_0.sigma+h_sigma, hat_A_G=gaussian_parameters_0.hat_A_G)
        # Get the starred values (i.e., stellar model's values at the position of the glitch)
        starred_values_h_sigma = glitch_wrt_model(gaussian_parameters_h_sigma.loc,unperturbed_stellar_model)
        # Create the Gaussian for the Gaussian parameters *with* increment
        N2_gaussian_h_sigma = get_new_N2(mask_gaussian=mask_gaussian,\
                                         stellar_model=unperturbed_stellar_model,\
                                         gaussian_modulation_buoyancy=starred_values_h_sigma['gaussian_modulation_buoyancy'],\
                                         gaussian_parameters=gaussian_parameters_h_sigma)
        # Obtain the Gaussian parameters w.r.t. the perturbed model
        gaussian_parameters_perturbed_model_h_sigma = gaussian_parameters_wrt_final_model(stellar_model=unperturbed_stellar_model,\
                                                                                          new_N2_gaussian=N2_gaussian_h_sigma,\
                                                                                          gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                          starred_values_wrt_original=starred_values_h_sigma,\
                                                                                          plot=False)

        # Compute the derivative of the Gaussian parameters w.r.t. the Gaussian's sigma 
        dhat_A_G_dsigma = ( gaussian_parameters_perturbed_model_h_sigma['hat_A_G'] - gaussian_parameters_perturbed_model_0['hat_A_G'] ) / h_sigma
        dsigma_dsigma   = ( gaussian_parameters_perturbed_model_h_sigma['sigma']   - gaussian_parameters_perturbed_model_0['sigma'] )   / h_sigma
        dloc_dsigma     = ( gaussian_parameters_perturbed_model_h_sigma['loc']     - gaussian_parameters_perturbed_model_0['loc']   )   / h_sigma

        # Collect the derivatives
        derivatives_wrt_sigma = { 'hat_A_G':dhat_A_G_dsigma, 'sigma':dsigma_dsigma, 'loc':dloc_dsigma }

        ### W.r.t gaussian hat A_G ### ---------------------------------------------------------

        gaussian_parameters_h_hat_A_G = Gaussian_parameters(loc=gaussian_parameters_0.loc, sigma=gaussian_parameters_0.sigma, hat_A_G=gaussian_parameters_0.hat_A_G+h_hat_A_G)
        # Get the starred values (i.e., stellar model's values at the position of the glitch)
        starred_values_h_hat_A_G = glitch_wrt_model(gaussian_parameters_h_hat_A_G.loc,unperturbed_stellar_model)
        # Create the Gaussian for the Gaussian parameters *with* increment
        N2_gaussian_h_hat_A_G = get_new_N2(mask_gaussian=mask_gaussian,\
                                           stellar_model=unperturbed_stellar_model,\
                                           gaussian_modulation_buoyancy=starred_values_h_hat_A_G['gaussian_modulation_buoyancy'],\
                                           gaussian_parameters=gaussian_parameters_h_hat_A_G)
        # Obtain the Gaussian parameters w.r.t. the perturbed model
        gaussian_parameters_perturbed_model_h_hat_A_G = gaussian_parameters_wrt_final_model(stellar_model=unperturbed_stellar_model,\
                                                                                            new_N2_gaussian=N2_gaussian_h_hat_A_G,\
                                                                                            gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                                                            starred_values_wrt_original=starred_values_h_hat_A_G,\
                                                                                            plot=False)

        # Compute the derivative of the Gaussian parameters w.r.t. the Gaussian's hat A_G 
        dhat_A_G_dhat_A_G = ( gaussian_parameters_perturbed_model_h_hat_A_G['hat_A_G'] - gaussian_parameters_perturbed_model_0['hat_A_G'] ) / h_hat_A_G
        dsigma_dhat_A_G   = ( gaussian_parameters_perturbed_model_h_hat_A_G['sigma']   - gaussian_parameters_perturbed_model_0['sigma'] )   / h_hat_A_G
        dloc_dhat_A_G     = ( gaussian_parameters_perturbed_model_h_hat_A_G['loc']     - gaussian_parameters_perturbed_model_0['loc']   )   / h_hat_A_G

        # Collect the derivatives
        derivatives_wrt_hat_A_G = { 'hat_A_G':dhat_A_G_dhat_A_G, 'sigma':dsigma_dhat_A_G, 'loc':dloc_dhat_A_G }

        # Collect all the derivatives
        derivatives = {'wrt_loc':derivatives_wrt_loc, 'wrt_sigma':derivatives_wrt_sigma, 'wrt_hat_A_G':derivatives_wrt_hat_A_G}

        # Return the derivatives
        return derivatives

    def button_clear_Gaussian(event):
        '''Clear the curves of the Gaussians from the plot'''
        clear_line2D(plots.figure, iplots['ax1']['gaussian_line2D'], plots.ax1, redraw=True)
        clear_line2D(plots.figure, iplots['ax1']['gaussian_modulation_line2D'], plots.ax1, redraw=True)
        clear_line2D(plots.figure, iplots['ax2']['gaussian_line2D'], plots.ax2, redraw=True)

    def overplot_eigenmode_wavelength_scale(event):
        '''Over plot estimations of mode's wavelength from a GYRE-generated file'''

        # Setup the root window.
        Tk().withdraw()
        # Show an dialog windows to select and return the path to the file.
        filename = askopenfilename()
        if Verbose:
            print('\tChosen file: ', filename)

        # Read the file
        header, data = tomso.gyre.load_mode(filename)

        # Find approximate positions where the change in sign occurs
        signs = np.sign(data['Rexi_r'])
        jumps = signs[1:] - signs[:-1]
        ind_nodes = np.where(np.abs(jumps) == 2)

        # Estimate nodes and wavelengths
        fractional_radius_for_nodes = data['x'][ind_nodes]
        fractional_radius_wavelength_estimation = fractional_radius_for_nodes[::2]

        # Restrict the wavelengths to the buoyancy cavity
        condition1 = stellar_model.fractional_r_buoyancy.min() < fractional_radius_wavelength_estimation
        condition2 = fractional_radius_wavelength_estimation < stellar_model.fractional_r_buoyancy.max()
        mask_fractional_radius_wavelength_estimation = np.logical_and(condition1, condition2)

        # Convert fractional radius to normalized buoyancy depth 
        normalized_buoyancy_depth_wavelength_estimation = np.interp(fractional_radius_wavelength_estimation, stellar_model.fractional_r_buoyancy, stellar_model.normalized_buoyancy_depth)

        # Plot as vertical lines
        iplots['ax1']['mode_line2D'], = axvlines(normalized_buoyancy_depth_wavelength_estimation, ax=plots.ax1, color='orange', label='Wavelength estimation', linewidth=0.5)
        iplots['ax2']['mode_line2D'], = axvlines(fractional_radius_wavelength_estimation[mask_fractional_radius_wavelength_estimation], ax=plots.ax2, color='orange', label='Wavelength estimation', linewidth=0.5)

        # Redraw
        plots.redraw()

    def clear_eigenmode_wavelength_scale(event):
        '''Clear the vertical lines for the wavelengths from the plot'''
        clear_line2D(plots.figure, iplots['ax1']['mode_line2D'], plots.ax1, redraw=True)
        clear_line2D(plots.figure, iplots['ax2']['mode_line2D'], plots.ax2, redraw=True)

    def get_the_mask_for_the_gaussian_in_ax1_interactively(stellar_model,\
                                                           min_normalized_buoyancy_depth,\
                                                           max_normalized_buoyancy_depth):
        '''Update the Gaussian mask'''

        clear_line2D(plots.figure, iplots['ax1']['gaussian_interval_line2D'], plots.ax1, redraw=False)
        clear_line2D(plots.figure, iplots['ax2']['gaussian_interval_line2D'], plots.ax2, redraw=False)
       
        # Print the selected interval
        print('\nSelected interval for the Gaussian (in units of total buoyancy depth):')
        print('\tmin = {:.3}, \tmax = {:.3}\n'.format(min_normalized_buoyancy_depth, max_normalized_buoyancy_depth))
        
        # Get the indices of the values within the selected span
        c1 = min_normalized_buoyancy_depth < stellar_model.normalized_buoyancy_depth
        c2 = stellar_model.normalized_buoyancy_depth < max_normalized_buoyancy_depth
        mask = np.logical_and(c1, c2)
        # Save the value ouside the scope of the function
        gaussian_interval['ax1']['mask'] = mask
        
        # Plot in red the selected span
        if min_normalized_buoyancy_depth != max_normalized_buoyancy_depth:
            iplots['ax1']['gaussian_interval_line2D'], = plots.ax1.plot(stellar_model.normalized_buoyancy_depth[mask], stellar_model.N2_buoyancy[mask], linestyle='None', marker='o', markerfacecolor='red', markeredgecolor='None', label='Range to be considered for the Gaussian')
            if stellar_model.type == 'adipls':
                iplots['ax2']['gaussian_interval_line2D'], = plots.ax2.plot(stellar_model.fractional_r_buoyancy[mask], stellar_model.A4_buoyancy[mask], linestyle='None', marker='o', markerfacecolor='red', markeredgecolor='None', label='Considered range')
            # Save the value ouside the scope of the function
            gaussian_interval['ax1']['ready'] = True
        else:
            # Save the value ouside the scope of the function
            gaussian_interval['ax1']['ready'] = False
        
        # Redraw
        plots.redraw()

    def plot_gaussian(mask_gaussian,\
                      normalized_buoyancy_depth,\
                      gaussian_modulation_buoyancy,\
                      modulation_scale,\
                      ax2x,\
                      ax2y):
        '''Plot the Gaussian in the two axes'''
        # Define some convenient temporal aliases
        modulation = gaussian_modulation_buoyancy[mask_gaussian]
        normalized_buoyancy_depth = normalized_buoyancy_depth[mask_gaussian]
        clear_line2D(plots.figure, iplots['ax1']['gaussian_line2D'], plots.ax1, redraw=False)
        clear_line2D(plots.figure, iplots['ax2']['gaussian_line2D'], plots.ax2, redraw=False)
        clear_line2D(plots.figure, iplots['ax1']['gaussian_modulation_line2D'], plots.ax1, redraw=False)
        # Get the y limits of the plot before ploting the modulation in order to put the same limits after plotting the modulation
        ylim = plots.ax1.get_ylim()
        # Plot the modulation in yellow
        iplots['ax1']['gaussian_modulation_line2D'], = plots.ax1.plot(normalized_buoyancy_depth, modulation_scale * modulation, linestyle='solid', marker='None', color='yellow', markerfacecolor='None', markeredgecolor='None', label=r'Modulation $\sqrt{N_0^\star/N_0\ r/r^\star}$ (arbitrary units)')
        # Set the y limits same as before plotting the modulation
        plots.ax1.set_ylim(ylim)
        # Plot the Gaussian
        iplots['ax1']['gaussian_line2D'], = plots.ax1.plot(normalized_buoyancy_depth, results['new_N2_gaussian'], linestyle='solid', marker='o', color='lime', markerfacecolor='lime', markeredgecolor='None', label=r'Modified $N^2$')        
        # Plot
        iplots['ax2']['gaussian_line2D'], = plots.ax2.plot(ax2x, ax2y, linestyle='solid', marker='o', color='lime', markerfacecolor='lime', markeredgecolor='None', label=r'Modified $A$')              
        # Redraw
        plots.redraw() 

    def increase_resolution_of_interval(min_normalized_buoyancy_depth,\
                                        max_normalized_buoyancy_depth,\
                                        extra_meshpoints,\
                                        overwrite=False):
        '''Increase the number of meshpoints in the model'''
        if stellar_model.type == 'adipls':
            # Find the correspond A0 (fractional radius)
            if np.all(np.diff(stellar_model.normalized_buoyancy_depth) > 0):
                fractional_r_max = np.interp(min_normalized_buoyancy_depth, stellar_model.normalized_buoyancy_depth, stellar_model.fractional_r_buoyancy)
                fractional_r_min = np.interp(max_normalized_buoyancy_depth, stellar_model.normalized_buoyancy_depth, stellar_model.fractional_r_buoyancy)
            else:
                fractional_r_max = np.interp(min_normalized_buoyancy_depth, stellar_model.normalized_buoyancy_depth[::-1], stellar_model.fractional_r_buoyancy[::-1])
                fractional_r_min = np.interp(max_normalized_buoyancy_depth, stellar_model.normalized_buoyancy_depth[::-1], stellar_model.fractional_r_buoyancy[::-1])

            stellar_model.increase_resolution(fractional_r_min, fractional_r_max, extra_meshpoints, overwrite=overwrite)

            if overwrite:
                plots.make_plots(stellar_model)
                plots.redraw()
            else:
                plots.make_plots(stellar_model.new_resolution)
                plots.redraw()

    def onselect(vmin, vmax):
        ''' Interactively select the intervals for the adding Gaussian '''

        if Verbose:
            comment = '[from the beginning of function onselect()]:\n| Prior values of the Gaussian parameters in units of normalized buoyancy before selecting an interval by mouse'
            gaussian_p.print(comment=comment)

        # Increase the resolution of an interval
        if (iinput.keystroke == 'ctrl+enter' and iinput.pressed_button == 1) \
        or (iinput.keystroke == 'enter' and iinput.pressed_button == 2):

            # To be on the safe side and force the Gaussian interval selection again
            gaussian_interval['ax1']['ready'] = False

            if Verbose:
                print('vmin =', vmin)
                print('vmax =', vmax)

            if vmin != vmax:
                # Setup the root window.
                Tk().withdraw()
                # Ask user for the number of new points 
                extra_meshpoints = simpledialog.askstring('Enter integer number', '\n     Number of new meshpoints to add     \n     between each original pair:     \n') or None
                # If the user does not input any value
                if extra_meshpoints is None:
                    pass
                else:
                    # If the user input a valid integer  
                    try:
                        extra_meshpoints = np.int(extra_meshpoints)
                        increase_resolution_of_interval(vmin, vmax, extra_meshpoints=extra_meshpoints)
                        iinput.new_resolution = True
                    except ValueError as e:
                        print('Harmless AttributeError :::', e)
                        print('| Input a valid integer.')
                        print('| Input value: {}'.format(extra_meshpoints))
            else:
                # Plot the original model
                iinput.new_resolution = False
                plots.make_plots(stellar_model)
                plots.redraw()

            # If the Gaussian interval had been previously selected, then get it again for the new resolution
            if gaussian_interval['ax1']['ready']:
                if iinput.new_resolution:
                    get_the_mask_for_the_gaussian_in_ax1_interactively(stellar_model.new_resolution,\
                                                                       min_normalized_buoyancy_depth=gaussian_interval['ax1']['min'],\
                                                                       max_normalized_buoyancy_depth=gaussian_interval['ax1']['max'])
                else:
                    get_the_mask_for_the_gaussian_in_ax1_interactively(stellar_model,\
                                                                       min_normalized_buoyancy_depth=gaussian_interval['ax1']['min'],\
                                                                       max_normalized_buoyancy_depth=gaussian_interval['ax1']['max'])

        # Select the interval for the gaussian by pressing the key 'enter' and then selecting it with the left button of the mouse
        if iinput.keystroke == 'enter' and iinput.pressed_button == 1:

            if Verbose:
                print('[from onselec()]:')
                print("|'if' block to select the Gaussian interval.")

            # Store the values in the nonlocal variables
            gaussian_interval['ax1']['min'] = vmin
            gaussian_interval['ax1']['max'] = vmax
            
            get_the_mask_for_the_gaussian_in_ax1_interactively(stellar_model.new_resolution if iinput.new_resolution else stellar_model,\
                                                               min_normalized_buoyancy_depth=gaussian_interval['ax1']['min'],\
                                                               max_normalized_buoyancy_depth=gaussian_interval['ax1']['max'])
        
        # Add the not normalized gaussian by pressing the key 'enter' and selecting the Gaussian parameters with the left button of the mouse
        # Alternatively, add the not normalized Gaussian by pressin the keys 'shift+enter' and selecting the Gaussian parameters with the rigth button of the mouse
        if (iinput.keystroke == 'shift+enter' and iinput.pressed_button == 1) \
        or (iinput.keystroke == 'enter' and iinput.pressed_button == 3):

            if vmin == vmax:
                clear_line2D(plots.figure, iplots['ax1']['gaussian_line2D'], plots.ax1)
                clear_line2D(plots.figure, iplots['ax2']['gaussian_line2D'], plots.ax2)
                clear_line2D(plots.figure, iplots['ax1']['gaussian_modulation_line2D'], plots.ax1)

            else:
                # Check if the interval where to add the Gaussian is already chosen
                if gaussian_interval['ax1']['ready']:

                    # Read the mouse-selected parameters

                    # (1) Gaussian location:
                    gaussian_p.update_loc(iinput.pressed_xdata)

                    # Get the starred values (i.e., stellar model's values at the position of the glitch)
                    starred_values = glitch_wrt_model(gaussian_p.loc,stellar_model.new_resolution if iinput.new_resolution else stellar_model)

                    # half FWHM -> (2) sigma:
                    gaussian_half_FWHM = vmax - vmin
                    gaussian_FWHM = 2 * gaussian_half_FWHM
                    gaussian_sigma = gaussian_FWHM / ( 2*np.sqrt( 2*np.log(2) ) )
                    gaussian_p.update_sigma(gaussian_sigma)

                    # new N2 -> (3) hat_A_G (keeping unchanged sigma):
                    results['gaussian_N2'] = iinput.released_ydata # Formerly called 'peak'
                    gaussian_amplitude = results['gaussian_N2']/starred_values['N2'] - 1
                    gaussian_hat_A_G = gaussian_amplitude * np.sqrt(2*np.pi) * gaussian_sigma
                    gaussian_p.update_hat_A_G(gaussian_hat_A_G)

                    results['new_N2_gaussian'] = get_new_N2(mask_gaussian=gaussian_interval['ax1']['mask'],\
                                                            stellar_model=stellar_model.new_resolution if iinput.new_resolution else stellar_model,\
                                                            gaussian_modulation_buoyancy=starred_values['gaussian_modulation_buoyancy'],\
                                                            gaussian_parameters=gaussian_p,\
                                                            comment='Output parameters for the mouse-selected Gaussian:')

                    ax2x, ax2y = get_data_for_ax2(mask_gaussian=gaussian_interval['ax1']['mask'],\
                                                  stellar_model=stellar_model.new_resolution if iinput.new_resolution else stellar_model)

                    # Do the plot
                    plot_gaussian(mask_gaussian=gaussian_interval['ax1']['mask'],\
                                  normalized_buoyancy_depth=stellar_model.new_resolution.normalized_buoyancy_depth if iinput.new_resolution else stellar_model.normalized_buoyancy_depth,\
                                  gaussian_modulation_buoyancy=starred_values['gaussian_modulation_buoyancy'],\
                                  modulation_scale=results['gaussian_N2']/starred_values['modulation_starred'],\
                                  ax2x=ax2x,\
                                  ax2y=ax2y)


                    results['ax2y'] = ax2y

                    _ = gaussian_parameters_wrt_final_model(stellar_model=stellar_model.new_resolution if iinput.new_resolution else stellar_model,\
                                                            new_N2_gaussian=results['new_N2_gaussian'],\
                                                            gaussian_mask=gaussian_interval['ax1']['mask'],\
                                                            starred_values_wrt_original=starred_values,\
                                                            plot=False)


    def window_info(_):
        '''Display information about how to use'''
        # Text separators
        sep1 = '--------------------------------------------------------------------------'
        sep2 = '================================================================'
        sep3 = ''
        # Info text
        title = 'Info'
        p0 = '>>                           Panels                           <<'
        p1 = ('> Left panel shows the buoyancy profile as parametrized in Eq.\n'
              '  (13) in Cunha et al. 2019. The Gaussian will be added to\n'
              '  this panel and then translated into the right panel')
        p2 = ('> Right panel shows the buoyancy profile as stored in the\n'
              '  AMDL model.')
        p3 = 'Both panels support the basic Matplotlib interactive tools'
        c0 = '>>                           Cursor                           <<'
        c1 = 'To enable recognition of the mouse buttons, first hit enter'
        c1 = '> The left button selects the interval to add the Gaussian.'
        c2 = ("> The right button defines the parameters of the Gaussian:\n\n"
              "    * The cursor's position when pressing defines the location.\n\n"
              "    * The cursor's position when releasing defines the peak.\n\n"
              "    * The distance between the pressing and releasing defines\n"
              "      the FWHM.")
        c3 = ('> The center button selects the interval where to increase the\n'
              '  number of meshpoints.')
        c4 = ("> The behaviour of the right button can be simulated by\n"
              "  the left button if this is preceded by the keystroke\n"
              " 'Shift+Enter'.\n\n"
              "  The behaviour of the center button can be simulated by the\n"
              "  left button if this is preceded by the keystroke 'Ctrl+Enter'.")
        b0 = '>>                           Buttons                          <<'
        b1 = ('> The button [Add Gaussian] adds a Gaussian using the values\n'
              '  enter in the text boxes. Notice that you have to hit Enter or\n'
              '  click on the text box to load the number that it displays.\n\n'
              '  If the button is pressed after adding a Gassian using the,\n'
              '  cursor then the loaded Gaussian parameters are the ones from\n'
              '  that Gaussian. ' )
        b2 = ('> The button [Add mode scale] plots estimations of the\n'
              '  wavelength from a GYRE mode file.')
        cb0 = '>>                          Check box                         <<'
        cb1 = ('> Mark the checkbox [w.r.t. final model] to make the button\n'
              '  [Add Gaussian] treat the loaded Gaussian parameters with\n'
              '  respect to the final model (with the added Gaussian) instead\n'
              '  of the original model.')
        e0 = '>>                           Exit                             <<'
        e1 = ('> When closing the plot, you will be asked in the terminal\n'
               '  whether to save the changes.')
        # Message to display    
        message = ('\n\n' + '\n\n'.join([sep2,p0,sep2,p1,sep3,p2,\
                                         sep2,c0,sep2,c1,sep3,c2,sep3,c3,sep3,c4,\
                                         sep2,b0,sep2,b1,sep3,b2,\
                                         sep2,cb0,sep2,cb1,\
                                         sep2,e0,sep2,e1]) + '\n\n')
        # Create blank windows
        window = tkinter.Tk()
        # Create title
        window.title(title)
        # Add content and organize its geometry
        tkinter.Label(window, text=message, justify = 'left').pack()
        # Modify canvas and organize its geometry
        tkinter.Canvas(window, width = 600, height=0).pack()
        # Display the window
        window.mainloop()

    def gaussian_parameters_wrt_final_model(stellar_model,\
                                            new_N2_gaussian,\
                                            gaussian_mask,\
                                            starred_values_wrt_original,\
                                            plot=True):

        # Copy the origianl
        new_N2_buoyancy = stellar_model.N2_buoyancy.copy()
        # Apply the modification
        new_N2_buoyancy[gaussian_mask] = new_N2_gaussian.copy()
        # Get new N
        new_N_buoyancy = np.sqrt(new_N2_buoyancy)
        # Get new omega_g
        new_omega_g_buoyancy = integrate.simps( stellar_model.L*new_N_buoyancy/stellar_model.r_buoyancy, stellar_model.r_buoyancy )
        # Buoyancy depth
        new_buoyancy_radius = integrate.cumtrapz(stellar_model.L*new_N_buoyancy/stellar_model.r_buoyancy, stellar_model.r_buoyancy, initial=0)
        # Take the outer turning point as the start or zero point
        new_buoyancy_depth = new_omega_g_buoyancy - new_buoyancy_radius
        # Normalize the buoyancy depth
        new_normalized_buoyancy_depth = new_buoyancy_depth/new_omega_g_buoyancy
        # Get new location
        new_gaussian_loc = np.interp(starred_values_wrt_original['r'], stellar_model.r_buoyancy, new_normalized_buoyancy_depth)
        
        # Generae the Gaussian-loke profile:
        scaled_gaussian_buoyancy = (new_N2_buoyancy/stellar_model.N2_buoyancy - 1) / starred_values_wrt_original['gaussian_modulation_buoyancy']
        
        # Functions to fit
        def fit_perturbated_model(x, hat_A_G, sigma):
            amplitude = hat_A_G / ( np.sqrt(2*np.pi) * sigma )
            return scaled_Gaussian(x=x, height=amplitude, loc=new_gaussian_loc, sigma=sigma )
        def fit_unperturbated_model(x, hat_A_G, sigma):
            amplitude = hat_A_G / ( np.sqrt(2*np.pi) * sigma )
            return scaled_Gaussian(x=x, height=amplitude, loc=gaussian_p.loc, sigma=sigma )

        # Fits
        hat_A_G_and_sigma_unperturbed_model, _ = curve_fit(fit_unperturbated_model, stellar_model.normalized_buoyancy_depth, scaled_gaussian_buoyancy, p0=np.array([gaussian_p.hat_A_G, gaussian_p.sigma]))
        try:
            hat_A_G_and_sigma_perturbed_model, _ = curve_fit(fit_perturbated_model, new_normalized_buoyancy_depth, scaled_gaussian_buoyancy, p0=np.array([gaussian_p.hat_A_G, gaussian_p.sigma]))
        except ValueError as e:
            print('ValueError :::', e)
            embed()

        gaussian_parameters_perturbed_model = {'hat_A_G':hat_A_G_and_sigma_perturbed_model[0], 'sigma':hat_A_G_and_sigma_perturbed_model[1], 'loc':new_gaussian_loc, 'omega_g_buoyancy':new_omega_g_buoyancy}
        gaussian_parameters_unperturbed_model = {'hat_A_G':hat_A_G_and_sigma_unperturbed_model[0], 'sigma':hat_A_G_and_sigma_unperturbed_model[1], 'loc':gaussian_p.loc, 'omega_g_buoyancy':stellar_model.omega_g_buoyancy}
        
        if plot:
            fig = plt.figure()
            ax = plt.axes()
            ax.plot(new_normalized_buoyancy_depth, scaled_gaussian_buoyancy, label='Perturbed model', color='red', linestyle='solid')
            ax.plot(new_normalized_buoyancy_depth, fit_perturbated_model(new_normalized_buoyancy_depth, *hat_A_G_and_sigma_perturbed_model), label=r'Perturbed model (estimated $\hat{A}_G$ and $\Delta_g$)', color='orange', linestyle='dashed')
            ax.plot(stellar_model.normalized_buoyancy_depth, scaled_gaussian_buoyancy, label='Unperturbed model', linestyle='solid', color='dodgerblue')
            ax.plot(stellar_model.normalized_buoyancy_depth, fit_unperturbated_model(stellar_model.normalized_buoyancy_depth, *hat_A_G_and_sigma_unperturbed_model), label='Unperturbed model (estimated $\hat{A}_G$ and $\Delta_g$)', linestyle='dashed', color='cyan')
            ax.axvline(new_gaussian_loc, color='lime', linestyle='dotted')
            ax.axvline(gaussian_p.loc, color='lime', linestyle='dotted')
            ax.invert_xaxis()
            plt.legend()
            plt.show()
        
        if Verbose:
            print('')
            print('[from gaussian_parameters_wrt_final_model():]')
            print('|')
            print('|- Original Gaussian parameters w.r.t. the unperturbed model')
            print('|Location = ', gaussian_p.loc)
            print('|Sigma = ', gaussian_p.sigma)
            print('|Hat A_G = ', gaussian_p.hat_A_G)
            print('|')
            print('|- Fit for the Gaussian parameters w.r.t. the unperturbated stellar model:')
            print('|Location (original parameter) = ', gaussian_parameters_unperturbed_model['loc'])
            print('|Sigma (fitted) = ', gaussian_parameters_unperturbed_model['sigma'])
            print('|Hat A_G (fitted) = ', gaussian_parameters_unperturbed_model['hat_A_G'])
            print('|omega_g  (original parameter) = ', gaussian_parameters_unperturbed_model['omega_g_buoyancy'])
            print('|')
            print('|- Fit for the Gaussian parameters w.r.t. the perturbated stellar model:')
            print('|Location (calculed) = ', gaussian_parameters_perturbed_model['loc'])
            print('|Sigma (fitted) = ', gaussian_parameters_perturbed_model['sigma'])
            print('|Hat A_G (fitted) = ', gaussian_parameters_perturbed_model['hat_A_G'])
            print('|omega_g  (calculed) = ', gaussian_parameters_perturbed_model['omega_g_buoyancy'])
            print('')

        # Return the Gaussian parameters hat A_G and sigma, and loc for the perturbed stellar model 
        return gaussian_parameters_perturbed_model

    def get_new_N2(mask_gaussian,\
                   stellar_model,\
                   gaussian_modulation_buoyancy,\
                   gaussian_parameters,\
                   comment='',\
                   Verbose=False):
        ''' Equeation (13) in the 2019 paper '''

        # Define some convenient temporal aliases
        N2 = stellar_model.N2_buoyancy[mask_gaussian] 
        modulation = gaussian_modulation_buoyancy[mask_gaussian]
        normalized_buoyancy_depth = stellar_model.normalized_buoyancy_depth[mask_gaussian]
        scaled_gaussian = scaled_Gaussian(x=normalized_buoyancy_depth, height=gaussian_parameters.amplitude, loc=gaussian_parameters.loc, sigma=gaussian_parameters.sigma )
        
        # Generate the new N2 profile
        new_N2_gaussian = N2 * ( 1 + modulation * scaled_gaussian )
        
        # Print the output values
        if Verbose:
            print('')
            print('[from get_new_N2():]')
            gaussian_p.print(comment=comment)
            print('')

        return new_N2_gaussian

    def get_data_for_ax2(mask_gaussian,\
                         stellar_model):

        # Prepare plot for the right axis
        if stellar_model.type == 'adipls':
            ax2x = stellar_model.fractional_r_buoyancy[mask_gaussian]
            ax2y = results['new_N2_gaussian'] * stellar_model.r_buoyancy[mask_gaussian] / stellar_model.g_buoyancy[mask_gaussian]
        elif stellar_model.type == 'gyre':
            ax2x = stellar_model.r_buoyancy[mask_gaussian]
            ax2y = results['new_N2_gaussian']

        return ax2x, ax2y

    # Gaussian parameters
    gaussian_p = Gaussian_parameters(loc=gaussian_loc_initial_text,\
                                     sigma=gaussian_sigma_initial_text,\
                                     hat_A_G=gaussian_hat_A_G_initial_text)

    # Results
    results = { 'gaussian_N2':     np.nan,\
                'new_N2_gaussian': np.nan,\
                'ax2y':            np.nan }

    # Text box inputs
    itextbox = { 'loc':     '',\
                 'sigma':   '',\
                 'hat_A_G': '',\
                 'N2peak':  '' }

    # To store the interval where to add the Gaussian
    gaussian_interval = { 'ax1':{ 'min':None, 'max':None, 'mask':None, 'ready':False },\
                          'ax2':{ 'min':None, 'max':None } }

    if interactive:

        # To store the changing plots
        iplots = { 'ax1':{ 'gaussian_line2D':            None,\
                           'gaussian_modulation_line2D': None,\
                           'gaussian_interval_line2D':   None,\
                           'mode_line2D':                None },\
                   'ax2':{ 'gaussian_line2D':            None,\
                           'gaussian_interval_line2D':   None,\
                           'mode_line2D':                None } }

        # To store interactive information by the user
        iinput = Interactive_variables()

        # ID connections to the plot visualization
        cids = { 'key':             plots.figure.canvas.mpl_connect('key_press_event', onPressKey),\
                 'button_pressed':  plots.figure.canvas.mpl_connect('button_press_event', onPressButton),\
                 'button_released': plots.figure.canvas.mpl_connect('button_release_event', onReleaseButton) }


        # Coordenates for the text boxes
        textboxes_coords = [  (0.23, 0.030, 0.05, 0.035), \
                              (0.37, 0.030, 0.05, 0.035), \
                              (0.49, 0.030, 0.05, 0.035), \
                              (0.60, 0.030, 0.05, 0.035)  ]

        # Particular parameters for each text box
        textbox_peak = Textbox(    function=read_peak,\
                                   prompt_text=r'(Peak)$_{\Delta_g}$ ($\mu$Hz$^2$)= ',\
                                   initial_text='{}'.format(gaussian_peak_initial_text),\
                                   coords=textboxes_coords[0] )
        
        textbox_hat_A_G = Textbox( function=read_hat_A_G,\
                                   prompt_text=r'$\hat{A}_G$ '+r'$(\omega_g)$ = ',\
                                   initial_text='{}'.format(gaussian_hat_A_G_initial_text), \
                                   coords=textboxes_coords[1] )
        
        textbox_sigma = Textbox(   function=read_sigma,\
                                   prompt_text=r'$\Delta_g$ '+r'$(\omega_g)$ = ',\
                                   initial_text='{}'.format(gaussian_sigma_initial_text),\
                                   coords=textboxes_coords[2] )
        
        textbox_loc = Textbox(     function=read_loc,\
                                   prompt_text=r'$\omega_g^{r \star} (\omega_g)$ = ',\
                                   initial_text='{}'.format(gaussian_loc_initial_text),\
                                   coords=textboxes_coords[3] )

        # Coordenates for the buttons
        buttons_coords = [  (0.740, 0.05, 0.1, 0.040), \
                            (0.740, 0.01, 0.1, 0.040), \
                            (0.845, 0.05, 0.1, 0.040), \
                            (0.845, 0.01, 0.1, 0.040),
                            (0.050, 0.01, 0.105, 0.040)  ]

        # Buttons
        button_add_mode_wavelenth = Button(     function=overplot_eigenmode_wavelength_scale,\
                                                text='Add mode scale',\
                                                coords=buttons_coords[0] )

        button_clear_mode_wavelength = Button(  function=clear_eigenmode_wavelength_scale,\
                                                 text='Clear mode scale',\
                                                coords=buttons_coords[1] )

        button_add_gaussian = Button(           function=add_gaussian,\
                                                text='Add Gaussian',\
                                                coords=buttons_coords[2] )

        button_clear_gaussian = Button(         function=button_clear_Gaussian,\
                                                text='Clear Gaussian',\
                                                coords=buttons_coords[3] )
    
        button_info = Button(         function=window_info,\
                                      text='Help',\
                                      coords=buttons_coords[4] )

      
        # Coordenates for the the check box
        checkbox_coords = [ (0.05, 0.05, 0.105, 0.04) ]

        # Checkbox
        checkbox_wrt_final = Checkbox(function=wrt_final_buoyancy_depth_check,\
                                      label='w.r.t. final model',\
                                      coords=checkbox_coords[0],\
                                      type=1)

        # Properties of the rectangle-span area-selector
        rect_props = dict(facecolor='cyan', alpha=0.20) 
        # Area selector
        span = mwidgets.SpanSelector(plots.ax1, onselect, 'horizontal',rectprops=rect_props, useblit=True) 

        # Make space for the text box below the plot
        plt.subplots_adjust(bottom=0.2)

        # Display plot
        plt.show()

        # Disconnect from the plot visualization
        for cid in cids.items():
            plots.figure.canvas.mpl_disconnect(cid)   

    # If not interactive
    else:

        # If increment in resolution in not interactive mode, first increase the resolution
        if bool(extra_meshpoints):
            increase_resolution_of_interval(gaussian_p.loc-4*gaussian_p.sigma,\
                                            gaussian_p.loc+4*gaussian_p.sigma,\
                                            extra_meshpoints=extra_meshpoints)
            mask_gaussian = np.zeros(stellar_model.new_resolution.normalized_buoyancy_depth.shape, dtype=bool)
        else:
            # Since it is the non-interactive mode, assume the interval to modify is the whole cavity
            mask_gaussian = np.zeros(stellar_model.normalized_buoyancy_depth.shape, dtype=bool)
        
        mask_gaussian[:] = True
        gaussian_interval['ax1']['mask'] = mask_gaussian
        gaussian_interval['ax1']['ready'] = True

        # Call the same function as the button "add gaussian". Emulate/simulate the "botton event" by inputing a None
        add_gaussian(None)

 
    if ( iinput.new_resolution if interactive else bool(extra_meshpoints) ):
        stellar_model.over_write_new_resolution()

    # Return the mask for the modified elements, the modified N2 and the modified variable to save into the stellar model
    return gaussian_interval['ax1']['mask'], results['new_N2_gaussian'], results['ax2y']

def modify_stellar_model(model,\
                         stellar_model_type,\
                         mask_buoyancy_cavity,\
                         new_data,\
                         save=True,\
                         filename_output=None):
    '''
    
    Purpose:

        Take the new values for the N2 profile and, consistently, modify the
        profile of the first adiabatid exponent (Gamma 1) as in the 2018 paper
        by Ball et. al. (https://arxiv.org/abs/1804.11153)

    Inputs:

        - model (string):

            Name of the stellar model file to be modified. The file can be
            either (1) an unformatted amdl file created by ADIPLS or (2) a
            MESA-GYRE formatted text file created by MESA.

        - stellar_model_type (string):

            Type of stellar model. Support for 'adipls' only.

        - mask_buoyancy_cavity (boolean array):

            Mask corresponding to the bouyancy cavity

        - new_data:

            output from the function "add_gaussian_interactively()"

    Optional inputs:

        - save (boolean):


    '''

    if stellar_model_type == 'adipls':
        # Parse the new data depending on the format of the stellar model
        mask_gaussian, _, A4_new = new_data
        if isinstance(model, Read_stellar_model):
            D = model.D
            A = model.A
        # Make a copy
        D_modified = D.copy()
        A_modified = A.copy()
        # Inverse of the pressure heitght scale from Hydrostatic equilibrium
        dlnP_dlnr = - ( A[:,3][mask_buoyancy_cavity][mask_gaussian] * A[:,2][mask_buoyancy_cavity][mask_gaussian] )
        # Inverse of the density height scale from the Brun-Vaisala frequency definition
        dlnrho_dlnr = - A[:,2][mask_buoyancy_cavity][mask_gaussian] - A[:,4][mask_buoyancy_cavity][mask_gaussian]
        # Mofify the adimensional N2. Make it equal to the fit
        A4_modified_temp_mask_buoyancy_cavity = A_modified[:,4][mask_buoyancy_cavity]
        A4_modified_temp_mask_buoyancy_cavity[mask_gaussian] = A4_new
        A_modified[:,4][mask_buoyancy_cavity] = A4_modified_temp_mask_buoyancy_cavity
        # Clear V_g from Gamma_1
        A_modified[:,2][mask_buoyancy_cavity][mask_gaussian] = A[:,2][mask_buoyancy_cavity][mask_gaussian] * A[:,3][mask_buoyancy_cavity][mask_gaussian]
        inv_Gamma_1 = (A_modified[:,4][mask_buoyancy_cavity][mask_gaussian] + dlnrho_dlnr) * dlnP_dlnr**-1
        # Set the new value for Gamma_1 to compensate for the change in N2
        A_modified[:,3][mask_buoyancy_cavity][mask_gaussian] = inv_Gamma_1**-1
        # Put back the (now modified) Gamma_1 in V_g
        A_modified[:,2][mask_buoyancy_cavity][mask_gaussian] = A_modified[:,2][mask_buoyancy_cavity][mask_gaussian] / A_modified[:,3][mask_buoyancy_cavity][mask_gaussian]
        # Save the modifications
        if save:
            if filename_output is None: 
                # Ask the user for a name
                filename_output = ask_output_name(model.filename+'.added_glitch')
            # If the user does not provide a name, then do not save the changes.
            if filename_output is None:
                print( '\t' + 'Model not saved.' )
            else:
                # Save the modified model 
                tomso.adipls.save_amdl(filename_output, D_modified, A_modified)
        else:
            print( '\t' + 'Model not saved.' )

def add_gaussian_A4(stellar_model,\
                    buoyancy_between,\
                    stellar_model_type,\
                    save=True,\
                    Verbose=False,\
                    G=6.67428e-08,\
                    gaussian_peak_initial_text='0.0',\
                    gaussian_hat_A_G_initial_text='0.0',\
                    gaussian_sigma_initial_text='0.0',\
                    gaussian_loc_initial_text='',\
                    interactive=True,\
                    filename_output=None,\
                    wrt_final_buoyancy_depth=False,\
                    extra_meshpoints=None): 
    '''

    Purpose:

        Add a Gaussian-like glitch to the N^2 profile of an existing stellar
        model in AMDL format, where N is the Brunt Vaisala frequency. The
        Gaussian-like function is parametrized as described in Cunha et al. 
        2019. The modification to the buoyancy frequency is done in a consistent
        way by a correspondent sole modification of the first adiabatic exponent
        as described in Ball et al. 2018, section 2.3.

        The proram offers an interactive GUI as well as a command line mode.

    Inputs:

        - stellar_model (string):

            String pointing to the AMDL file to modify.

        - buoyancy_between (float,float):

            Upper and lower limits where to search for the bouyancy cavity.
            Values in unit of fractional radius. For instance, if the star has
            a convective core, a sensible choice is (0.0,0.5).

            It is necessary to define this range in order to not mix the
            buoyancy cavity with other regions where the buoyancy frequency is
            also different from 0, for example, near-surface regions.

        - stellar_model_type (string):

            So far, there is only support for stellar model generated by ADIPLS.
            Other types may be incorporated in the future.

    Optional inputs:

        - save=True (boolean):

            Whether or not to save the changes at the end.

        - Verbose=False (boolean):

            Print additional information on the terminal. Often useful for
            debugging.

        - G=6.67428e-08 (float):

            Gravitational constant defined as positive in units of dyne cm**2 g**-1.
            Value taken from https://tomso.readthedocs.io/en/latest/adipls.html

        - gaussian_peak_initial_text='' (float):

            Value to show as initial text in the text box in the GUI. It 
            corresponds to the peak of the added Gaussian in units of Hz^2.

            Only allowed in the if interactive=True

        - gaussian_hat_A_G_initial_text=0.00343776 (float):

            Value to show as initial text in the text box in the GUI. It
            corresponds to the amplitud of the added Gaussian in units of total
            buoyancy depth.

            If the program is used in the non-interactive version, this is the
            value used as the amplitude of the Gaussian

        - gaussian_sigma_initial_text=0.00048 (float):

            Value to show as initial text in the text box in the GUI. It
            corresponds to the sigma of the added Gaussian in units of total
            buoyancy depth.

            If the program is used in the non-interactive version, this is the
            value used as the sigma of the Gaussian

        - gaussian_loc_initial_text=0.00568 (float):

            Value to show as initial text in the text box in the GUI. It
            corresponds to the position of the added Gaussian in units of total
            buoyancy depth.

            If the program is used in the non-interactive version, this is the
            value used as the amplitude of the Gaussian

        - interactive=True (boolean):

            Whether to go into the interactive GUI.

        - filename_output=None (string):

            Name to save the modified AMDL stellar model

        - wrt_final_buoyancy_depth=False (boolean):

            Whether to specify the parameters of the Gaussian with respect to
            the unperturbed model (the initial one) or the perturbed model (the
            final one with the added Gaussian).

        - extra_meshpoints=False (integer):

            Only valid when interactive=False. Number of new meshpoints to place
            in between the original meshpoints of the model. This only takes
            effect about the center of the Gaussian (3 sigma in units of total
            buoyancy depth). 

    Outputs:

        The program will output a new AMDL stellar model

    Dependences:

        Make sure you have installed the following Python packages:

            - Numpy (https://www.numpy.org/)
            - Matplotlib (https://matplotlib.org/)
            - Scipy (https://www.scipy.org/)
            - TOMSO (https://tomso.readthedocs.io/en/latest/#)  

    Example:

        Let 'model.amdl' be an AMDL stellar model generated by ADIPLS.

        There are 2 following ways to use this program.

        1) As an executable program from the terminal. For instance:

                $ ./add_gaussian_A4 model.amdl -t adipls -b 0.0,0.5 -l 0.00568 -w 0.00048 -A 0.00343776 -i yes -s yes -v yes

           When used as executable, the inputs are given through the
           following flags (-t and -b are mandatory):

                -t: stellar_model_type (mandatory)

                -b: buoyancy_between (mandatory)

                -v: Verbose
                
                -s: save
                
                -G: Gravitational constans
                
                -p: gaussian_peak_initial_text

                -l: gaussian_loc_initial_text

                -w: gaussian_sigma_initial_text

                -A: gaussian_hat_A_G_initial_text

                -i: interactive

                -o: filename_output

                -f: wrt_final_buoyancy_depth

                -r: extra_meshpoints

                -h: This will show this help.

           If you are using a different python kernel from the default one
           (e.g. conda), then it may be better to specify the interpreter as,
           for instance:

                $ python3 add_gaussian_A4 model.amdl -t adipls -b 0.0,0.5

        2) Within Python by importing the function. For instance:
            
                >>> from add_gaussian_A4 import add_gaussian_A4
                >>> add_gaussian_A4('model.amdl', (0.0,0.5), 'adipls', gaussian_hat_A_G_initial_text='0.00343776', gaussian_loc_initial_text='0.00568', gaussian_sigma_initial_text='0.00048')

    '''

    # Read the stellar model
    model = Read_stellar_model(filename=stellar_model,\
                               created_by=stellar_model_type,\
                               buoyancy_cavity_between=buoyancy_between,\
                               G=G,\
                               verbose=Verbose)

    # Make a plot of N2 vs buoyancy depth and a second plot vs the fractional radius
    plots = Make_plots(model)

    # Modify N2
    new_model = add_gaussian_interactively(model,\
                                           plots,\
                                           gaussian_peak_initial_text=gaussian_peak_initial_text,\
                                           gaussian_hat_A_G_initial_text=gaussian_hat_A_G_initial_text,\
                                           gaussian_sigma_initial_text=gaussian_sigma_initial_text,\
                                           gaussian_loc_initial_text=gaussian_loc_initial_text,\
                                           Verbose=Verbose,\
                                           interactive=interactive,\
                                           wrt_final_buoyancy_depth=wrt_final_buoyancy_depth,\
                                           extra_meshpoints=extra_meshpoints)

    # Record the modified N2 profile and compensate by modifying the Gamma1 as in Ball et. at. 2018
    if save:
        modify_stellar_model(model,\
                             stellar_model_type,\
                             model.mask_buoyancy,\
                             new_model,\
                             filename_output=filename_output)

# If called as an script from the terminal
if __name__ == "__main__":
    import sys
    import re
    # Optional inputs conveniently displayed.
    optional_input = { 'gaussian_peak_initial_text':'',\
                       'gaussian_hat_A_G_initial_text':0.001750,\
                       'gaussian_sigma_initial_text':0.000263,\
                       'gaussian_loc_initial_text':0.050000,\
                       'save':False,\
                       'Verbose':True,\
                       'interactive':True,\
                       'filename_output':None,\
                       'wrt_final_buoyancy_depth':False,\
                       'extra_meshpoints':None,\
                       'G':6.67428e-08,\
                       'interactive':True,\
                       'filename_output':None}
    # Parse the user's command-line input
    filename, modeltype, buoyancy_between, flags = parse_user_command_line_input(sys.argv)
    # Over write the optional inputs in 'flags' with the optional inputs in 'kwargs'
    optional_input = { **optional_input, **flags }
    add_gaussian_A4(filename, buoyancy_between, modeltype, **optional_input)