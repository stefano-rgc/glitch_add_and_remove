import pandas as pd  
import numpy as np    
import tomso.adipls
from scipy import integrate
from copy import deepcopy
from scipy.interpolate import interp1d

class Read_stellar_model:
    ''' Collect and organize the information about the stellar model'''

    def __init__(self,\
                 filename,\
                 created_by,\
                 buoyancy_cavity_between,\
                 notes='',\
                 G=6.67428e-08,\
                 l=1,\
                 verbose=False):

        '''
        Purpose:

            Read the relevant data from a stellar model created by ADIPLS
            (an unformated AMDL file).

        Inputs:

            - filename (string):

                Name pointing to the stellar model file to be modified.

            - created_by (string):

                Type of stellar model. So far, only support for 'adipls'.

            - buoyancy_cavity_between (float, float):

                Tuple indicating the upper and lower limits of the interval
                where to search for the buoyancy cavity (in fractional radius).

        Optional Inputs:

            - notes='' (string):

                Comments.

            - G=6.67428e-08 (float):

                Gravitational constant defined as positive in units of dyne cm**2 g**-1.
                Value taken from https://tomso.readthedocs.io/en/latest/adipls.html

            - l=1 (integer):

                Angular degree mode. It enters the calculation of the total size
                of the buoyancy cavity.

            - verbose=False:

                Print extra information useful for debugging.
        '''

        # Support only for stellar model in ADIPLS format
        if created_by != 'adipls':
            raise ValueError('Please make sure to set the variable "created_by" to "adipls". Support for MESA-GYRE stellar model may be incorporated in the future.')
 
        # Basic info:
        self.type = created_by
        self.filename = filename
        self.G = G
        self.notes = notes
        self.l = l
        self.L = np.sqrt( self.l * ( self.l + 1 ) )
        self.buoyancy_cavity_between = buoyancy_cavity_between

        self.read_data_form_file(verbose=verbose)
        self.create_profiles(verbose=verbose)

    def add_comment(self, comment, verbose=False):
        self.notes = self.notes + '\n\n' + comment
        if verbose:
            print( "The following comment has been added as a new line to the attribute 'notes':" )
            print( "|New comment: {}\n".format(comment) )

    def show_notes(self): # This can be improved by decorators
        print(self.notes)

    def read_data_form_file(self, verbose=False):
        if self.type == 'adipls':
            D, A = tomso.adipls.load_amdl(self.filename)
            self.D = D
            self.A = A
            if verbose:
                print('AMDL attributes "A" and "D" have been set.\n')

    def set_adipls_new_data(self, D, A, verbose=False):
        if self.type != 'adipls':
            raise ValueError("The attribute 'type' is not 'adipls'.\n")
        self.D = D
        self.A = A
        if verbose:
            print( 'AMDL attributes "A" and "D" have been overwritten.')
            print( 'Profiles will be overwritten accordingly.\n' )
        self.create_profiles(verbose=verbose)

    def create_profiles(self, verbose=False):
        if self.type == 'adipls':
            self.A0 = self.A[:,0]
            self.A1 = self.A[:,1]
            self.A2 = self.A[:,2]
            self.A3 = self.A[:,3]
            self.A4 = self.A[:,4]
            self.A5 = self.A[:,5]
            # Derivated quantities
            R, r, g = tomso.adipls.amdl_get(['R', 'r', 'g'], self.D, self.A, G=self.G)
            self.r = r # Units: cm
            self.R = R # Units: cm
            self.g = g # Units: cm/s**2
            self.fractional_r = self.r/self.R
            # Note this does not include the turning points but only the region in between
            condition1 = self.A4 > 0
            condition2 = self.fractional_r > self.buoyancy_cavity_between[0]
            condition3 = self.fractional_r < self.buoyancy_cavity_between[1]
            ind_buoyancy_cavity = np.where( (condition1) & (condition2) & (condition3) )
            self.ind_buoyancy = ind_buoyancy_cavity
            mask_buoyancy_cavity = np.zeros(self.A4.shape, dtype=bool)
            mask_buoyancy_cavity[ind_buoyancy_cavity] = True
            self.mask_buoyancy = mask_buoyancy_cavity
            # For convenience, also define the above variables only for the buoyancy cavity only
            self.A0_buoyancy = self.A0[self.mask_buoyancy]
            self.A1_buoyancy = self.A1[self.mask_buoyancy]
            self.A2_buoyancy = self.A2[self.mask_buoyancy]
            self.A3_buoyancy = self.A3[self.mask_buoyancy]
            self.A4_buoyancy = self.A4[self.mask_buoyancy]
            self.A5_buoyancy = self.A5[self.mask_buoyancy]
            self.r_buoyancy = self.r[self.mask_buoyancy]
            self.g_buoyancy = self.g[self.mask_buoyancy]
            self.fractional_r_buoyancy = self.fractional_r[self.mask_buoyancy]
            # Get the dimensional N2 and N frequencies within the buoyancy cavity
            # For consistency, keep the sufix "_buoyancy"
            self.N2_buoyancy = self.g_buoyancy * self.A4_buoyancy / self.r_buoyancy #Uunits: Hz**2
            self.N_buoyancy = np.sqrt(self.N2_buoyancy) # Units Hz
            # Total buoyancy depth. Note that the integral do not include the immediate vicinity of the turning points
            self.omega_g_buoyancy = integrate.simps( self.L*self.N_buoyancy/self.r_buoyancy, self.r_buoyancy )
            # Buoyancy depth
            self.buoyancy_radius = integrate.cumtrapz(self.L*self.N_buoyancy/self.r_buoyancy, self.r_buoyancy, initial=0)
            # Take the outer turning point as the start or zero point
            self.buoyancy_depth = self.omega_g_buoyancy - self.buoyancy_radius
            # Normalize the buoyancy depth
            self.normalized_buoyancy_depth = self.buoyancy_depth/self.omega_g_buoyancy
            # Normalize the buoyancy radius
            self.normalized_buoyancy_radius = self.buoyancy_radius/self.omega_g_buoyancy
            if verbose:
                # Print notice
                print( "The following atributes has been set:" )
                print( "| 'A0', 'A1', 'A2', 'A3', 'A4', 'A5'" )
                print( "| 'R', 'r', 'g', 'fractional_r'" )
                print( "| 'ind_buoyancy', 'mask_buoyancy'" )
                print( "| 'A0_buoyancy', A1_buoyancy', 'A2_buoyancy'" )
                print( "| 'A3_buoyancy', 'A4_buoyancy', A5_buoyancy'" )
                print( "| 'r_buoyancy', 'g_buoyancy', 'fractional_r_buoyancy'" )
                print( "| 'N2_buoyancy', 'N_buoyancy', 'omega_g_buoyancy'" )
                print( "| 'buoyancy_radius', 'buoyancy_depth'")
                print( "| 'normalized_buoyancy_depth''and''normalized_buoyancy_radius'")

    def increase_resolution(self, fractional_r_min, fractional_r_max, extra_meshpoints, verbose=False, overwrite=False):
        try: # To avoid any copy into itself more than one. 
            del self.new_resolution
        except AttributeError as e:
            if verbose:
                print('Harmless AttributeError:::', e)
        finally: # Make a copy of the model to later change its mesh
            self.new_resolution = deepcopy(self)

        self.new_resolution.extra_meshpoints = extra_meshpoints

        if self.type == 'adipls':
            # Mask indicating where to increase the resolution
            c1 = self.fractional_r >= fractional_r_min
            c2 = self.fractional_r <= fractional_r_max
            ind_higher_resolution = np.where( c1 & c2 )
            mask_higher_resolution = np.zeros(self.A0.shape, dtype=bool)
            mask_higher_resolution[ind_higher_resolution] = True
            # Fraction distances between one element and the next one in A0 (fractional radius)
            fractions = np.arange( 0.0, 1.0, 1/(1+extra_meshpoints) )
            # The fraction 0.0 corresponds to the original point
            fractions = fractions[1:]
            # Get the differences between fractional radius
            A0_diff_higher_resolution = np.diff(self.A0[mask_higher_resolution])
            # Initialize variable where to save the new data
            A0_higher_resolution = np.array([])
            for fraction in fractions:
                A0_higher_resolution = np.concatenate( [ A0_higher_resolution, self.A0[mask_higher_resolution][:-1] + A0_diff_higher_resolution * fraction ] )
            # Add the original points
            A0_higher_resolution = np.concatenate( [ A0_higher_resolution, self.A0 ] )
            A0_higher_resolution = np.sort( A0_higher_resolution )
            # Get the other variables in higher resolution
            A_given_A0 = interp1d(self.A0, self.A, axis=0)
            A_higher_resolution = A_given_A0(A0_higher_resolution)
           
        if overwrite:
            self.set_adipls_new_data(self.D,A_higher_resolution, verbose=verbose)
        else:
            self.new_resolution.ind = ind_higher_resolution
            self.new_resolution.mask = mask_higher_resolution
            self.new_resolution.set_adipls_new_data(self.D,A_higher_resolution, verbose=verbose)

    def over_write_new_resolution(self, verbose=False):
        self.set_adipls_new_data(self.D, self.new_resolution.A, verbose=verbose)
        del self.new_resolution
        print('The resolution of the model has changed.\n')

    def save(self, output_name, verbose=False):
        if self.type == 'adipls':
            tomso.adipls.save_amdl(output_name, self.D, self.A)
            if verbose:
                print("Stellar model saved under the name '{}'\n".format(output_name) )

def read_famdl(filename):
    '''
    Radius is famdl[:,0], Gamma1 is famdl[:,3], etc
    (as in section FAMDL from https://www.astro.up.pt/corot/ntools/docs/CoRoT_ESTA_Files.pdf)
    '''
    famdl = pd.read_fwf(filename, skiprows=2, widths=[20,20,20,20])    
    famdl = famdl.values.reshape(-1)                                                                                                              
    famdl = famdl[~np.isnan(famdl)].reshape(-1,6)
    return famdl
