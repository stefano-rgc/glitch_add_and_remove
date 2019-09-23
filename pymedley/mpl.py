import numpy as np
import matplotlib.pyplot as plt

################################################################################

def clear_line2D(Figure, Line2D, Axes, redraw=True):
    '''Check if the given Line2D is part of the axis. If yes, remove it from the axis'''
    if np.isin(Line2D, Axes.lines):
        Axes.lines.remove(Line2D)            
        if redraw == True:
            Axes.legend(loc='best', ncol=1, framealpha=0.5, fontsize=10)
            Figure.canvas.draw()

################################################################################

def axvlines(xs, ax=None, **plot_kwargs):
    '''

    Purpose:

        Plot vertical lines as one single matplotlib object Line2D.

    Source:

        https://stackoverflow.com/questions/24988448/how-to-draw-vertical-lines-on-a-given-plot-in-matplotlib

    Author: 

        Peter (Stackoverflow user)

    Original comments by Peter:

        Draw vertical lines on plot
        :param xs: A scalar, list, or 1D array of horizontal offsets
        :param ax: The axis (or none to use gca)
        :param plot_kwargs: Keyword arguments to be passed to plot
        :return: The plot object corresponding to the lines.

    '''

    if ax is None:
        ax = plt.gca()

    xs = np.array( (xs, ) if np.isscalar(xs) else xs, copy=False )
    lims = ax.get_ylim()
    x_points = np.repeat( xs[:, None], repeats=3, axis=1 ).flatten()
    y_points = np.repeat( np.array( lims + (np.nan, ) )[None, :], repeats=len(xs), axis=0 ).flatten()
    
    plot = ax.plot(x_points, y_points, scaley = False, **plot_kwargs)
    
    return plot
