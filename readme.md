# Bastr - A parallel solver for 2D/3D Burgers turbulence

This code is adapted from [astr](https://github.com/lcs27/astr) code - it keeps its parallel framework. However, the calculation process is of no relation with the astr code: we use a pseudo-spectral method.

## Input - Run
The input file should contail two files:
1. An initial flow field `datin/flowini2d.h5` or `datin/flowini3d.h5`, generated from pp(postprocess) module of astr or Bastr, containing flow field of u1, u2 and u3(for 3D) 

2. An `datin/input.dat` file with following format:
```
########################################################################
#                     input file of BASTR code                          #
########################################################################

# ia,ja,ka: The size of grid. All non-zero--3D mode; ka=0--the 2d mode 
2048,2048,0

# ref_t,reynolds: Reference variables for viscosity
273.15d0, 5000.d0

# maxstep,deltat,lwsequ,feqwsequ,lwspectra,feqwspe, kmax: output configurations
2000, 1.d-4, t, 2000, t, 200, 200

# forcemethod, targetenergy, lproject: Forcing and calculation
2,0.075d0,t
```
List of parameters:
- `maxstep` : Maximum running step
- `deltat` : $\Delta t$ of each step
- `lwsequ` : t/f to open sequential output
- `feqwsequ` : frequency(nstep) of sequential output
- `lwspectra` : t/f to open spectral output
- `feqwspe` : frequency(nstep) of spectral output
- `kmax` : maximum k for spectral output
- `forcemethod` : 0 = free decaying, 1 = linear forcing, 2 = forcing on a band
- `targetenergy` : valable only if `forcemethod` on 
- `lproject` : only project to dilatational part


## Dependencies
This code requires a dependency on
- mpi ([mpich](https://www.mpich.org/) or [openmpi](https://www.open-mpi.org/) is ok!)
- [hdf5](https://www.hdfgroup.org/download-hdf5/)
- [fftw](https://www.fftw.org/)

## Installation procedure
If your computer is already equipped with ASTR code? - Feel free! Just run make and you can compile the Bastr code.

If you begin from zero:
*TODO*