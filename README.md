# Physical context

The asymptotic analysis of the pulsation equations predicts regularly-spaced periods for gravity
modes. These periods get deflected when a sharp variation in the buoyancy frequency, associated with
strong chemical gradients, is present. These features are known as glitches and their study is one of
the few methods allowing us to probe the physical conditions in a localized region inside a star.
Theoretical works (e.g., Cunha et al. 2015, 2019) provide the analytical bases to interpret the
signatures of glitches of arbitrary amplitude in seismic observables. However, testing the analytical
predictions is not straightforward, because finding a model with a glitch with specific properties (e.g.,
position, amplitude, and width) implies, in the best-case scenario, to explore the space of stellar
parameters until the desired glitch is attained. This approach has clear disadvantages like (1) no
warranty to generate the pursued glitch, (2) time consuming, and (3) preclusion of a direct comparison
of the effect of different glitches on the seismic data.

This is a Python program that allows to artificially add/remove a Gaussian-like
glitch into/from the buoyancy profile of an already existing stellar model, that way solving points (1-
3). Comparison between the sharpness of the glitch with respect to the local wavelength of gravity
modes is also provided as well as the possibility to remesh the model to allow for a very narrow glitch.
The modification to the buoyancy frequency is done in a consistent way by a correspondent sole
modification of the first adiabatic exponent as described in Ball et al. 2018, section 2.3. The program
counts with an interactive GUI as well as a command line mode.

# Remove a glitch from a stellar model:

Current version can operate only with stellar models stored as unformatted AMDL
files, which is the format that the Aarhus adiabatic pulsation code (ADIPLS, Chistensen-
Dalsgaard 2008) uses.

Given an unformated AMDL stellar model, the program modifies the variable A4 (square of the dimensionless buoyancy frequency) within it
by fitting a polynomial function. The code is intended to remove sharp variations in the A4 profile and therefore
get rid of buoyancy glitches. The code has a GUI which allows the user to explore the A4 profile and
interactively select the interval to fit a polynomial while excluding a
smaller interval within itself (the latter is suppoused to contain the
glitch). The modification in the A4 profile is done in a consistent way by a
correspondent sole modification of the first adiabatic exponent, as
described in the section 2.3 of Ball et al. 2018 ("Surface effects on
the red giant branch").



The plots are generated by the Python package Matplotlib that provides interactive tools like zooming, dragging,
and rescaling, which are all useful for data exploration. The figure below shows and example of the removal of a glitch

![alt text](https://github.com/stefano-rgc/glitch_add_and_remove/blob/master/exemplary_images/remove_glitch.gif)

Short explanation:

- First, it zooms in on the glitch-like feature.
- Second, it selects the interval to be considered for the fit (red curve).
- Third, it selects the interval to be ecluded from the fit (black vertical lines).
- Fourth, it performs the polynomial fit by selecting a degree and then clicking the button (green curve).

After closing the windows, the program will prompt the user whether to save the changes and a name to output the new AMDL model. The following summary will be display on the console:

```
=========================================
Interval for the fit:
xmin = 0.0143, 	xmax_fit = 0.0313
=========================================
Interval excluded from the fit:
xmin = 0.0192, 	xmax = 0.0252
=========================================
Polynomial degree of the fit = 3
=========================================
```

To run the example, try to excute the following line in the console:

```
$python3 polypatch_A4.py  model.amdl 
```


# Add a Gaussian-like glitch to a stellar model:
