# Brief Explanation of Main Modules.

## reach.py
reach.py is responsible for constructing the initial set from a Model object and computing the resulting flowpipe. The computation returns a FlowPipe object which is passed to a FlowPipePlotter object for plotting.

## bundle.py
bundle.py contains all of the routines required to transform a bundle. This includes constructing the proper polynomials relevant to finding the offsets for the directions matrix and finding the maximum and minimum Bernstein coefficients for those polynomials.

## flowpipe.py
Contains the plotting routines interfacing with matplotlib and scipy. The routines are responsible for plotting the projections and phase plots of the reachable sets.

## bernstein.py
Responsible for the conversion of real polynomials into their counterparts expressed in the Bernstein basis. Furthermore, it extracts out the maximum and minimum Bernstein coefficients from the converted polynomial.

## parallelotope.py
Defines the Parallelotope class which contains all routines required to compute properties of the relevant parallelotope.

## model.py
Defines the parent Model class.

## lputil.py
Linear programming utilites using swiglpk.
