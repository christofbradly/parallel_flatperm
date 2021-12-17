# parallel_flatperm
Parallel implementation of flatPERM algorithm for lattice polymer simulations.

This program simulates self-avoiding walk models on a regular lattice using the flatPERM algorithm[^1] with parallelised implementation[^2].
Valid lattices are square, simple cubic, hypercubic up to 10 dimensions, triangular and hexagonal.
Most functionality is for square, simple cubic and triangular, limited funcitonality for others.
The algorithm can simulate self-avoiding walks, self-avoiding trails, neighbour-avoiding walks or grooves, where valid.
Available parameters for the flatPERM histogram include length n (required), near-neighbours m, stiffness s, contacts with surface a, height above surface h, multiply-visited sites (trails only) c and d.
Some others may be valid depending on the particular model chosen.
The algorithm can also be run in a fixed weight configuration where a parameter's conjugate Boltzmann weight is provided as initial input.
The weights of each sample are updated using this input weight in addition to the Rosenbluth & Rosenbluth weight. 
In this configuration the associated parameter is typically not flattened over in the histogram, but can be forced to do so.
The parallelisation uses a shared memory paradigm, with multiple threads simulating local chains and updating a shared global histogram. Race conditions are ignored with no significant loss of efficiency (cf. Wang-Landau algorithm[^3]).

## Usage
Options detailing the model to be simulated are input as command line arguments. Type
```main --help```
for further details.
Scripts for running on clusters are not provided.

Additional requirements:
 - GCC gfortran: Tested with v8.3. Support for at least Fortran 2003 standard is required.
 - [HDF5 v1.10+](https://www.hdfgroup.org/downloads/hdf5/): (earlier versions are not compatible with extended precision floats). Use provided h5fc compiler.
 - [mt_Stream](http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html): Fortran implementation of MT19937 random number generator for multiple parallel streams. Module and object files are included here but this package usually requires recompilation from source.
 - [fhash](https://github.com/jl2922/fhash): Fortran implementation of a generic hash table. I have included a version with some necessary modifications. In particular, the original node removal functinoality does not work.


## Notes on variables
The main program has a set of data arrays (histograms) and a loop over a defined number of tours.
Within each tour there is a loop that attempts to grow a walk step-by-step up to some maximum length `Nmax`, using the flatPERM algorithm to choose/commit/prune/enrich valid steps.
The state of the walk is stored in a number of variables:
- `chain`: an integer array of the occupied lattice sites in order
- `visits`: a hash table for querying which steps are occupying a lattice point
- `micro`: a short vector storing the microcanonical parameters of the state, used to index the histograms
- `history_`: a variety of arrays recording the value of several state variables in the order that the algorithm proceeds through state space (not the SAW)
- `i_hist`: an index for the history arrays

The use of history arrays is in lieu of a recursive implementation.
The files `microcanonical.f90` and `radius.f90` contain wrappers for calculating microcanonical parameters and metric quantities.
`savedata.f90` handles the HDF5 routines for saving the histograms.
Other files are for initialising the program, handling input arguments etc.




[^1]: [Prellberg & Krawczyk, *PRL*, 2004](https://doi.org/10.1103/PhysRevLett.92.120602)
[^2]: [Campbell & Janse van Rensburg, *JPA*, 2020](https://doi.org/10.1088%2F1751-8121%2Fab8ff7)
[^3]: [Zhan, *Comput. Phys. Commun.*, 2008](http://www.sciencedirect.com/science/article/pii/S0010465508001392)

