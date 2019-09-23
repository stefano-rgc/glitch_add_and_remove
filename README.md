# glitch_add_and_remove

## Context

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

## Program to remove a glitch:

![alt text]()

## Program to add a Gaussian-like glitch:
