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

The current version can operate only with stellar models stored as unformatted AMDL
files, which is the format that the Aarhus adiabatic pulsation code (ADIPLS, Chistensen-
Dalsgaard 2008) uses.

Given an unformated AMDL stellar model, the program modifies the variable A4 (square of the dimensionless buoyancy frequency) within it
by fitting a polynomial function. The code is intended to remove sharp variations in the A4 profile and therefore
get rid of buoyancy glitches. The code has a GUI which allows the user to explore the A4 profile and
interactively select the interval to fit a polynomial while excluding a
smaller interval within itself (the latter is supposed to contain the
glitch). The modification in the A4 profile is done in a consistent way by a
correspondent sole modification of the first adiabatic exponent, as
described in section 2.3 of Ball et al. 2018 ("Surface effects on
the red giant branch").

The code forces the fit to pass through the two last and two
first data point in the selected interval. This, in order to achieve a 
smoother matching between the data and the fit. As a result, **the degree 
of the fit must be 3 or higher**. This is done by means of Lagrange multipliers.

The plots are generated by the Python package Matplotlib that provides interactive tools like zooming, dragging,
and rescaling, which are all useful for data exploration. The figure below shows an example of the removal of a glitch

![alt text](https://github.com/stefano-rgc/glitch_add_and_remove/blob/master/exemplary_images/remove_glitch.gif)

Short explanation:

- First, it zooms in on the glitch-like feature.
- Second, it selects the interval to be considered for the fit (red curve).
- Third, it selects the interval to be included from the fit (black vertical lines).
- Fourth, it performs the polynomial fit by selecting a degree and then clicking the button (green curve).

After closing the windows, the program will prompt the user whether to save the changes and a name to output the new AMDL model. The following summary will be displayed on the console:

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

### Run the program

Make sure you are using Python3 and have installed the following Python packages:

  - Numpy (https://www.numpy.org/)
  - Matplotlib (https://matplotlib.org/)
  - Scipy (https://www.scipy.org/)
  - TOMSO (https://tomso.readthedocs.io/en/latest/#)  


To reproduce the example, execute the following line in the console:

```
$python3 polypatch_A4.py  model.amdl 
```
Extra documentation can be accessed from within Python via

```
>>> help(polypatch_A4)
```

# Add a Gaussian-like glitch to a stellar model:

> **The parametrization of the Gaussian in this code follows an earlier version of Cunha et al. 2019. Such parametrization, eq. (13) within the paper, was updated in its latest version (https://arxiv.org/abs/1909.04966) and the current version of this code still needs to be accordingly. As temporary mitigation, such update in the parametrization does not invalidate this code.**

The current version can operate only with stellar models stored as unformatted AMDL
files, which is the format that the Aarhus adiabatic pulsation code (ADIPLS, Chistensen-
Dalsgaard 2008) uses.

Given an unformatted AMDL stellar model, the program modifies the variable A4 (square of the dimensionless buoyancy frequency) within it by adding a Gaussian-like glitch accordingly to Eq. (13) in Cunha et al. 2019. The modification in the A4 profile is done in a consistent way by a correspondent sole modification of the first adiabatic exponent, as described in section 2.3 of Ball et al. 2018 ("Surface effects on the red giant branch"). The input parameters are the amplitude, width and position of the Gaussian. 

The interface of the program consist of two panels:

- The right panel shows the buoyancy profile as stored in the AMDL model.
- The left panel shows the buoyancy profile as parametrized in Eq. (13) in Margarida et al. 2019. **The Gaussian will be added to this panel and then translated into the right panel**

Plots are generated by the Python package Matplotlib that provides interactive tools like zooming, dragging,
and rescaling, which are all useful for data exploration.

The program provides the following three functionalities:

### (1) Addition of the glitch

![alt text](https://github.com/stefano-rgc/glitch_add_and_remove/blob/master/exemplary_images/add_glitch1.gif)

Short explanation:

- First, it zooms in on the region where the glitch will be added.
- Second, it selects the interval where the glitch will be added. (red curve).
- Third, it uses the cursor to add and define the parameters of the Gaussian (green curve).
- Fourth, Check how the changes look in the AMDL file (right panel).

### (2) Comparison with the mode wavelength

![alt text](https://github.com/stefano-rgc/glitch_add_and_remove/blob/master/exemplary_images/add_glitch2.gif)

Short explanation:

- First, it loads a **file of modes generated by GYRE** by clicking the button.
- Second, it changes the horizontal scale from logarithmic to linear to show how the wavelength is constant in the left panel (distances between orange vertical lines) 
- Third, it zooms in on the region where the glitch will be added.
- Fourth, it adds the glitch in the same way as section (1).


### (3) Remesh of the model (narrow glitches)

![alt text](https://github.com/stefano-rgc/glitch_add_and_remove/blob/master/exemplary_images/add_glitch3.gif)

Short explanation:

- First, it zooms in until the resolution of the meshpoints is visible.
- Second, adds a glitch which is undersampled.
- Third, it selects the interval where to remesh the model. Then a dialog pops up asking how many new meshpoints add between each original pair.
- Fourth, it selects the interval where the glitch will be added. Then adds it not by the cursor but by clicking the button (**done this way because, if not new Gaussian parameters are specificated, then the program uses the ones from the previous Gaussian, in this case, the undersampled Gaussian**).

After closing the windows, the program will prompt the user whether to save the changes and a name to output the new AMDL model.

### Run the program

Make sure you are using Python3 and have installed the following Python packages:

  - Numpy (https://www.numpy.org/)
  - Matplotlib (https://matplotlib.org/)
  - Scipy (https://www.scipy.org/)
  - TOMSO (https://tomso.readthedocs.io/en/latest/#)

To reproduce the examples, execute the following line in the console:

```
$python3 -W ignore add_gaussian_A4.py model.amdl -t adipls -b 0,0.5 
```

The option **'-W ignore'** is optional and hides the warning prints.

The program offers multiple options. A more complete example would be the following;

```
$python3 -W ignore add_gaussian_A4.py model.amdl -t adipls -l 0.00568 -w 0.00048 -A 0.00343776 -i yes -b 0,0.5 -s yes -v not -o model.amdl.modified -f not -i yes
```

**In particular, the flag '-i' (whose possible values are yes/not) determines if the program runs interactively or not**. For the meaning of the other flags, see the documentation of the program via 

```
$python3 add_gaussian_A4.py -h
```
