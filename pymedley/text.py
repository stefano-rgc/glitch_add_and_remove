import numpy as np

################################################################################

def parse_text_for_number_or_blank(text,\
                                   variable,\
                                   variable_name='Variable',\
                                   Verbose=True, ):
    '''
    Purpose:
        
        Recive a 'text', if it is a blank (''), then t

    Inputs:

        - text (string):

            String that will be check to be (1) a number, (2) a blank ('') or
            (3) something else. If (1), then the string is converted to a numpy
            float and assigned to the variable 'variable'. If (2), the a string
            blank ('') is assigned to the variable 'variable'. If (3), no new
            value is assigned to the variable 'variable'.

        - Variable (any type):

            

        - variable_name (string):

            String that will be printed as a message to the user to refer to the
            variable 'variable'


    Outputs:


        - If (1), 'numpy.float(text)'. If (2), 'text'. If (3), 'variable'.

    '''

    if text != '':
        try:
            if Verbose:
                print('{} = {}'.format(variable_name, variable))
            return np.float(text)
        except ValueError:
            print('\tInput a valid number for the {}.'.format(variable_name))
            return variable
    else:
        if Verbose:
            print('{} = {}'.format(variable_name, variable))
        return text

################################################################################

def ask_output_name(tentative_name):
    '''Ask the user for a name while showing a suggestion'''
    new_name = input('Output name: [{}]. Type n and Enter to cancel.\n'.format(tentative_name))
    if new_name == 'n': # In case the user cancels
        return None
    else:
        if new_name == '': # In case the user doesn't enter a different name
            return tentative_name
        else:
            return new_name

################################################################################

relate_flags_and_variable_names = {'-t':'stellar_model_type',\
                                   '-b':'buoyancy_cavity_between',\
                                   '-s':'not_save',\
                                   '-v':'Verbose',\
                                   '-G':'G',\
                                   '-p':'gaussian_peak_initial_text',\
                                   '-A':'gaussian_hat_A_G_initial_text',\
                                   '-w':'gaussian_half_FWHM_initial_text',\
                                   '-l':'gaussian_location_initial_text'}

valid_values_for_the_flags = {'-t':{'adipls', 'amdl', 'gyre'},\
                              '-s':{'not', 'yes', 'true', 'false'},\
                              '-v':{'not', 'yes', 'true', 'false'}}

convert_flag_values = {'yes':True,\
                       'true':True,\
                       'not':False,\
                       'false':False}

def parse_user_command_line_input(number_of_mandatory_inputs,\
                                  command_line_text_list,\
                                  relate_flags_and_variable_names,\
                                  valid_values_for_the_flags,\
                                  convert_flag_values,\
                                  case_fold_for_flag_values=True):
    '''

    Purpose:

    Inputs

        command_line_text_list

            Example:

                'input1 input2 -a flag1 -b flag2 -c flag3'

                command_line_text_list = [ 'input1', 'input2', '-a', 'flag1', '-b', 'flag2', '-c', 'flag3' ]

        number_of_mandatory_inputs

    '''

    valid_command_line = True

    # Split the command line into mandatory arguments and flags
    mandatory_inputs = command_line_text_list[number_of_mandatory_inputs:]
    flags = args[number_of_mandatory_inputs::2]
    flags_value = args[number_of_mandatory_inputs+1::2]

    # Check the number of mandatory inputs match the command line input text
    difference = number_of_mandatory_inputs - len(mandatory_inputs)
    if difference != 0:
        print( '**It seems {} mandatory inputs are missing**'.format(difference) )
        return not valid_command_line

    # Handle the case
    if case_fold_for_flag_values:
        flags_value = [ str(fv).casefold() for fv in flags_value ]

    # Check if the mandatory inputs are not flags
    for s in mandatory_inputs:
        if len(s) == 2 and s[0] == '-':
            print( '**The mandatory input {} looks like a flag.**'.format(s) )
            return not valid_command_line
    
    # Number of arguments
    narg = len(command_line_text_list)

    # Initialize flags
    valid_input_line = True


    # If an odd number of command-line inputs is given
    if odd_number(narg):

        # Check if help is requested
        if narg == 3 and command_line_text_list[2] == '-h':

            


        else:

            print("**Wrong number of command line inputs. Please see the help.**")
            return None, None, None, not valid_input_line
    
    # If an even number of command-lines inputs is given
    else:
        # Assume the flags start form the third comand-line argument
        flags = command_line_text_list[2::2]
        flags_value = command_line_text_list[3::2]
        # Handle the case
        flags_value = [ str(fv).casefold() for fv in flags_value ]
        
        # Check the name of the flags are a subset of the avaiblable flags
        if set(flags).issubset( relate_flags_and_variable_names.keys() ):
            
            # Check for repeated flags
            if len( set(flags) ) < len(flags):
                print("**Repeated flags.**")
                return None, None, None, not valid_input_line
            
            # If no repeated flags
            else:
                # Check for valid flag values
                for f,fv in zip(flags, flags_value):
            
                    # Some flags like '-G' need this check
                    if f in valid_values_for_the_flags:

                        if not fv in valid_values_for_the_flags[f]:
                            print( "**Please, select from among '{}' a valid input for the flag '{}'.**".format( valid_values_for_the_flags[f], f ) )
                            return None, None, None, not valid_input_line
            
                # If all the flags ara valid
                flags_and_values = { f:v for f,v in zip(flags,flags_value) }
                # Get the filename to read
                filename = str(command_line_text_list[1])

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
                        print("**Please, specify the type of stellar model (amdl or gyre) by the flag '-t'.**")
                        return None, None, None, not valid_input_line 
            
                # Check if extrama of the interval (in units of fractional radius) where to search for the buoyancy cavity are set
                if '-b' in flags:
                    # Extract the flag from the dictionary
                    temp = flags_and_values.pop('-b')
                    # Check that the input are two numbers separated by a comma
                    temp = temp.split(',')
                    if len(temp) == 2:
                        try:
                            buoyancy_cavity_between = ( np.float(temp[0]), np.float(temp[1]) )
                        except ValueError:
                            print("**Input two valid numbers in units of fractional radius separated by a comma and no space to bound the buoyancy cavity.**")
                            return  None, None, None, not valid_input_line
                    else:
                        print("**Input two valid numbers in units of fractional radius separated by a comma and no space to bound the buoyancy cavity.**")
                        return  None, None, None, not valid_input_line
                else:
                    print("**The flag '-b' is a mandatory flag. Please see help.**")
                    return  None, None, None, not valid_input_line


                # Check if a value for the gravitational constant is set
                if '-G' in flags:
                    try:
                        flags_and_values['-G'] = np.float( flags_and_values['-G'] )
                    except ValueError:
                        print("**Input a valid number for the gravitational constant in cgs units (e.g. 6.67428e-08)**")
                        return  None, None, None, not valid_input_line

                # Check if the location for the Gaussian-like function is set
                if '-l' in flags:
                    try:
                        flags_and_values['-l'] = np.float( flags_and_values['-l'] )
                    except ValueError:
                        print("**Input a valid number for the location of the Gaussian-like function in units of normalized buyancy depth**")
                        return  None, None, None, not valid_input_line    

                # Check if the value of the peak of the Gaussian-like function is set
                if '-p' in flags:
                    try:
                        flags_and_values['-p'] = np.float( flags_and_values['-p'] )
                    except ValueError:
                        print("**Input a valid number for the value of the peak of the Gaussian-like function in units of Hz^2**")
                        return  None, None, None, not valid_input_line                    

                # Check if the value of the half of the FWHM of the Gaussian-like function is set
                if '-w' in flags:
                    try:
                        flags_and_values['-w'] = np.float( flags_and_values['-w'] )
                    except ValueError:
                        print("**Input a valid number for the half of the FWHM of the Gaussian-like function in units of normalized buoyancy depth**")
                        return  None, None, None, not valid_input_line                    

                # Check if the parameter for the amplitud of the Gaussian-like function is set (as described in Margarida et. al.'s 2019 paper)
                if '-A' in flags:
                    try:
                        flags_and_values['-A'] = np.float( flags_and_values['-A'] )
                    except ValueError:
                        print("**Input a valid number for the amplitud of the Gaussian-like function in units of normalized buoyancy depth (as described in Margarida et. al.'s 2019 paper)**")
                        return  None, None, None, not valid_input_line                    

                # Convert the user command-line flag-value input into proper Python-language variables
                for f in ['-s', '-v']:
                    if f in flags_and_values:
                        flags_and_values[f] = convert_flag_values[ flags_and_values[f] ]

                # Convert the user command-line flag input accordingly to relate_flags_and_variable_names
                flags_and_values = { relate_flags_and_variable_names[f]:fv for f,fv in flags_and_values.items() }
                
                # Return parsed command line
                return filename, stellar_model_type, buoyancy_cavity_between, flags_and_values, valid_input_line                    
        
        # If unknown flag set
        else:

            print("**Wrong flag as input. Please see the help.**")
            return None, None, None, not valid_input_line
