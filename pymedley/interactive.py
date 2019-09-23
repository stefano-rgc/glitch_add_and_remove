import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import tkinter
from pymedley.mpl import clear_line2D
from pymedley.math import polyfit_with_fixed_points

### Classes

class Textbox:
    '''Create a text entry on a Matplotlib figure'''
    def __init__(self, function, prompt_text, initial_text, coords, label_pad=0):
        self.function = function
        self.prompt_text = prompt_text
        self.initial_text = initial_text
        self.coords = coords
        self.axis = plt.axes(self.coords)
        self.textbox = mwidgets.TextBox(self.axis, self.prompt_text, initial=self.initial_text, label_pad=label_pad) 
        self.textbox.on_submit(self.function)

class Button:
    '''Create a button on a Matplotlib figure'''
    def __init__(self, function, text, coords):
        self.function = function
        self.text = text
        self.coords = coords
        self.axis = plt.axes(self.coords)
        self.button = mwidgets.Button(self.axis, self.text)    
        self.button.on_clicked(self.function)

class Checkbox:
    '''Create a check box on a Matplotlib figure'''
    def __init__(self, function, label, coords, type=1, color=None):
        self.function = function
        self.label = label
        self.coords = coords
        self.axis = plt.axes(self.coords)
        self.checkbox = mwidgets.CheckButtons(self.axis, [self.label])    
        self.checkbox.on_clicked(self.function)
        if type == 1: # Check box this a cross
            self.checkbox.rectangles[0].set_y(0)
            self.checkbox.rectangles[0].set_x(0)
            self.checkbox.rectangles[0].set_height(1)   
            self.checkbox.rectangles[0].set_width(0.2)   
            self.checkbox.lines[0][0].set_xdata([0.0,0.2])
            self.checkbox.lines[0][0].set_ydata([0.0,1.0])
            self.checkbox.lines[0][1].set_xdata([0.0,0.2])
            self.checkbox.lines[0][1].set_ydata([1.0,0.0])
        if type == 2: # Check box with color only
            if color is None:
                color = 'lime' 
            self.checkbox.rectangles[0].set_facecolor(color)
            self.checkbox.rectangles[0].set_fill(False)
            self.checkbox.lines[0][0].set_linewidth(0)
            self.checkbox.lines[0][1].set_linewidth(0)
            self.checkbox.rectangles[0].set_y(0)
            self.checkbox.rectangles[0].set_x(0)
            self.checkbox.rectangles[0].set_height(1)   
            self.checkbox.rectangles[0].set_width(0.2)   

### Functions

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

def fit_polynomial(figure, axes, line2D, verbose=False):
    '''
    
    Description:

        Given a plot of Y vs X, select an interval of X and fit a polynomial to
        the data.

        Optionally, select a sub-interval of X to be ignored by the fit.

        IMPORTANT: The code forces the fit to pass through the two last and two
        first data point in the selected interval. This in order to achieve a 
        smoother mathing between the data and the fit. As a result, the degree 
        of the fit must be 3 or higher. This is done by means of Lagrange
        multipliers.

        The figure given as input can contain a legend, title, labels, etc.

    Inputs:

        figure:

            It refers to the figure object from Matplotlib.

        axes: 

            It refers to the axes object from Matplotlib.

        line2D:

            It refers to the line2D object from Matplotlib.

    Optional inputs:

        verbose:

            Print additional information on the terminal. Often useful for
            debugging.

    Output:

        The program returns a tupe where the first element is a mask which 
        contains the indices of the orifinal X corresponding to the fitted 
        interval. The second element are the new values of that interval.  

    Example:

        import numpy as np
        import matplotlib.pyplot as plt
        from pymedley.interactive import fit_polynomial

        # Sample data
        x = np.arange(-2,2,0.05)
        y = (x-1)*(x-0.5)*(x-0.2)*(x-0.1)*(x-2)*(x-1.5)*(x+1)*(x+1.3)*(x+1.1)*(x+2)

        # Making a plot
        figure = plt.figure()
        axes = plt.axes()
        line2d, = plt.plot(x,y) # Mind the comma!

        # Fitting a polynomial interactively
        mask, fit = fit_polynomial(figure, axes, line2d)

        # Integration of the fit values into the original data
        y[mask] = fit

    '''

    def window_info(_):
        '''Display information about how to use'''
        # Text separators
        l1 = '--------------------------------------------------------------------------'
        l2 = '=========================================================================='
        # Info text
        title = 'Info'
        t0 = '>> To enable recognition of the left and right buttons, first hit enter <<'
        t1 = '> The left button selects the interval for the fit by dragging on the plot.'
        t2 = '> The right button selects an interval excluded to by excluded from the fit.'
        t3 = '> Press the button [Add fit] to generate the fit.'
        t4 = '> Mark the checkbox [Not refresh fit] to overplot different fits.'
        t5 = '> When closing plot, you will be asked whether to save the changes.'
        # Message to display    
        message = ('\n\n' + '\n\n'.join([l2,t0,l2,t1,l1,t2,l1,t3,l1,t4,l1,t5,l1]) + '\n\n')
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

    def window_error():
        '''Display information about how to use'''
        title = 'Error'
        text = 'Please input an integer equal or greater than 3.'
        message = ('\n\n' + text + '\n\n')
        window = tkinter.Tk()
        window.title(title)
        tkinter.Label(window, text=message, justify = 'center').pack()
        tkinter.Canvas(window, width = 600, height=0).pack()
        window.mainloop()

    def read_keystroke(event):
        '''Get the pressed key over the axes during plot visualization'''
        ivar['keystroke'] = event.key
        if verbose:
            print( "ivar['keystroke'] = {}".format(ivar['keystroke']) )

    def read_button(event):
        '''Get the pressed button over the axes during plot visualization'''
        ivar['pressed_button'] = event.button
        if verbose:
            print("ivar['pressed_button']", ivar['pressed_button'])

    def read_polynomial_degree(text):
        '''Get the input text in the textbox during plot visualization'''
        if text.isdigit() and np.int(text)>2 :
            ivar['polynomial_degree'] = np.int(text)
        else:
            window_error()
            # Keep the previous value and display it on the check box.
            interface['textbox_polynomial_degree'].textbox.set_val(ivar['polynomial_degree'])
        if verbose:
            print( "ivar['polynomial_degree'] = {}".format(ivar['polynomial_degree']) )

    def button_add_fit(event):
        '''Add the fit to the plot'''
        nonlocal fit_function
        
        if fit_interval_ready:
            # If a previous fit, clear it
            clear_line2D(figure, lines2D['fit'], axes, redraw=False)
            if not ivar['not_refresh_flag']:
                clear_line2D(figure, lines2D['fit_denser'], axes, redraw=False)

            # Get current colors in the axis to not repeat them
            colors = [ l.get_color() for l in axes.get_lines() ]
            
            # Fit
            if mask['nonfit'] is not None:
                if interval['xmin_nonfit'] > interval['xmin_fit'] and interval['xmax_nonfit'] < interval['xmax_fit']:
                    mask_fit_minus_nonfit = np.logical_xor(mask['fit'], mask['nonfit'])
                else:
                    print('**The non fit interval is outside the fit interval. It will be ignored.**')
                    mask_fit_minus_nonfit = mask['fit'].copy() 
            else:
                mask_fit_minus_nonfit = mask['fit'].copy() 

            xfit = x[mask_fit_minus_nonfit]
            yfit = y[mask_fit_minus_nonfit] 

            # Force the fit to match the fist two data points and the last two data points to promote a smoother match with the original
            xfit_fixed = np.concatenate([ [xfit[0]], [xfit[1]], [xfit[-2]], [xfit[-1]] ])
            yfit_fixed = np.concatenate([ [yfit[0]], [yfit[1]], [yfit[-2]], [yfit[-1]] ])
            # Fit
            fit_function = np.poly1d( polyfit_with_fixed_points(ivar['polynomial_degree'], xfit, yfit, xfit_fixed, yfit_fixed) )
            color = 'lime' if not 'lime' in colors else None
            lines2D['fit'], = axes.plot(x[mask['fit']], fit_function(x[mask['fit']]), linestyle='None', marker='o', color=color, markerfacecolor=color, markeredgecolor='None')
            
            # Plot a denser x to evidence possible wiggles between the original data
            points_in_between = np.arange(0.0, 1.0, 0.1)
            x_mask_diff = np.diff(x[mask['fit']])
            x_mask_denser = np.array([])
            for shift in points_in_between:
                x_mask_denser = np.concatenate( [ x_mask_denser, x[mask['fit']][:-1] + x_mask_diff * shift ] )
            # Add the last point
            x_mask_denser = np.concatenate( [ x_mask_denser, np.array( [ x[mask['fit']][-1] ] ) ] )
            x_mask_denser = np.sort( x_mask_denser )
            color = lines2D['fit'].get_color()
            lines2D['fit_denser'], = axes.plot(x_mask_denser, fit_function(x_mask_denser), linestyle='solid', marker=None, color=color, label='Polynomial fit n={}'.format(ivar['polynomial_degree']))
            list_fit_denser_line2D.append(lines2D['fit_denser'])

        # Redraw
        axes.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)
        figure.canvas.draw()

    def button_clear_fit(event):
        '''Clear the fits curves from the plot'''
        nonlocal list_fit_denser_line2D
        nonlocal fit_function
        clear_line2D(figure, lines2D['fit'], axes, redraw=True)
        for line2D in list_fit_denser_line2D:
            clear_line2D(figure, line2D, axes, redraw=True)
        fit_function = None

    def check_box_switch(label):
        '''Invert the flag of not refresh'''
        ivar['not_refresh_flag'] = not ivar['not_refresh_flag']
        interface['checkbox'].checkbox.rectangles[0].set_fill(ivar['not_refresh_flag'])
        figure.canvas.draw()

    def onselect(vmin, vmax):
        '''Select the interval for the fit and the interval to be excluded in the fit'''

        nonlocal fit_interval_ready

        # Interval to fit
        # Activate by pressing Enter and then using the left button
        if ivar['keystroke'] == 'enter' and ivar['pressed_button'] == 1:
            clear_line2D(figure, lines2D['fit_interval'], axes, redraw=False)
            # Store the values
            interval['xmin_fit'] = vmin
            interval['xmax_fit'] = vmax
            # Print the interval
            print('Interval for the fit:')
            print('xmin = {:.3},\t xmax = {:.3}\n'.format(vmin, vmax))
            # Get the indices of the values within the selected span
            condition1_fit = interval['xmin_fit'] < x 
            condition2_fit = x < interval['xmax_fit']
            mask['fit'] = np.logical_and(condition1_fit, condition2_fit)
            # Plot in red the selected span as an aditional Line2D object in lines
            if interval['xmin_fit'] != interval['xmax_fit']:
                lines2D['fit_interval'], = axes.plot(x[mask['fit']], y[mask['fit']], linestyle='None', marker='o', markerfacecolor='red', markeredgecolor='None', label='Fit interval')
                fit_interval_ready = True
            else:
                fit_interval_ready = False

        # Interval to exclude
        # Activate by pressing Enter and then using the right button OR pressing Shift+Enter and the using left button
        if (ivar['keystroke']=='shift+enter' and ivar['pressed_button']==1) or (ivar['keystroke']=='enter' and ivar['pressed_button']==3):
            clear_line2D(figure, lines2D['v1'], axes, redraw=False)
            clear_line2D(figure, lines2D['v2'], axes, redraw=False)
            interval['xmin_nonfit'] = vmin
            interval['xmax_nonfit'] = vmax
            # Print the interval
            print('Interval to be excluded in the fit:')
            print('xmin = {:.3},\t xmax = {:.3}\n'.format(vmin, vmax))
            # Get the indices of the values within the selected span
            condition1_nonfit = interval['xmin_nonfit'] < x 
            condition2_nonfit = x < interval['xmax_nonfit']
            mask['nonfit'] = np.logical_and(condition1_nonfit, condition2_nonfit)
            # Plot in black the selected span
            if interval['xmin_nonfit'] != interval['xmax_nonfit']:
                lines2D['v1'] = axes.axvline(interval['xmin_nonfit'], label='Nonfit interval', linestyle='dashed', color='black')
                lines2D['v2'] = axes.axvline(interval['xmax_nonfit'], linestyle='dashed', color='black')

        # Redraw
        axes.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)
        figure.canvas.draw()

    ### Initialize variables

    # Curves in the plot
    lines2D = {'fit':None,
               'fit_interval':None,
               'fit_denser':None,
               'v1':None,
               'v2':None}
    # Masks
    mask = {'fit':None,
            'nonfit':None}
    # Variable to store the fit
    fit_function = None
    # Interval's extrema
    interval = {'xmin_fit':None,
                'xmax_fit':None,
                'xmin_nonfit':None,
                'xmax_nonfit':None}
    # Interactive variables
    ivar = {'keystroke':None,
            'pressed_button':None,
            'polynomial_degree':3,
            'not_refresh_flag':False}
    # List where to store curves of different fits
    list_fit_denser_line2D = list()
    # Flags
    fit_interval_ready = False

    # Make space for the interface of buttons
    plt.subplots_adjust(bottom=0.2)

    # Get the data
    x = line2D.get_xdata()
    y = line2D.get_ydata()

    # Connect ID to the plot visualization
    cid_key = figure.canvas.mpl_connect('key_press_event', read_keystroke)
    cid_button = figure.canvas.mpl_connect('button_press_event', read_button)

    # Biuld interface: Buttons, checkbox and textbox
    height = 0.04
    width  = 0.1
    position = [(0.65, 0.07)]
    button_addfit             = Button(  function=button_add_fit,\
                                         text='Add fit',\
                                         coords=[0.65, 0.07, width, height] )
    button_clearfit           = Button(  function=button_clear_fit,\
                                         text='Clear fit',\
                                         coords=[0.76, 0.07, width, height] )
    button_help               = Button(  function=window_info,\
                                         text='Help',\
                                         coords=[0.76, 0.02, width, height] )
    checkbox                  = Checkbox(function=check_box_switch,\
                                         label='Not refresh fit',\
                                         coords=[0.65, 0.02, width, 0.040],\
                                         type=2)
    textbox_polynomial_degree = Textbox( function=read_polynomial_degree,\
                                         prompt_text='Polynomial degree for the fit = ',\
                                         initial_text='{}'.format(ivar['polynomial_degree']),\
                                         coords=[0.45, 0.05, 0.05, 0.05])

    interface = {'button_addfit':button_addfit,
                 'button_clearfit':button_clearfit,
                 'button_help':button_help,
                 'checkbox':checkbox,
                 'textbox_polynomial_degree':textbox_polynomial_degree}


    # Properties of the rectangle-span area-selector
    rect_props = dict(facecolor='cyan', alpha=0.20) 
    # Area selector
    span = mwidgets.SpanSelector(axes, onselect, 'horizontal',rectprops=rect_props) 

    # Display plot in a maximazed window
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    plt.show()

    # Disconnect from the plot visuzlization
    for cid in [ cid_key, cid_button ]:
        figure.canvas.mpl_disconnect(cid)

    # Print the selected intervals
    print('\n')
    print('=========================================')
    print('Interval for the fit:')
    print('xmin = {:.3}, \txmax = {:.3}'.format(interval['xmin_fit'], interval['xmax_fit']))
    print('=========================================')
    print('Interval excluded from the fit:')
    print('xmin = {:.3}, \txmax = {:.3}'.format(interval['xmin_nonfit'], interval['xmax_nonfit']))
    print('=========================================')
    print('Polynomial degree of the fit = {}'.format(ivar['polynomial_degree']))
    print('=========================================')
    print('\n')

    # Return velues
    if fit_function != None:
        return mask['fit'], fit_function(x[mask['fit']])
    else:
        return None, None