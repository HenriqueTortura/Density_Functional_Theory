# Density Functional Theory in Simple Systems

## Build Instructions
Compile main module for Fortran (in my case, f95 calls for GNU Fortran)
```bash
f95 -c dft.f95
```
and the module for Python (currently, I could only make it work with two separated files)
```bash
f2py -c -m pydft pydft.f95
```
## Running
Call for hydrogen and/or helium solutions (open files to edit parameters)
```bash
f95 -o HydrogenAtom HydrogenAtom.f95 dft.o
f95 -o HeliumAtom HeliumAtom.f95 dft.o
./HydrogenAtom 
./HeliumAtom
```
Sweep hydrogen eigenvalues
```bash
f95 -o Sweep Sweep.f95 dft.o
./Sweep
```
Make hydrogen plots with `Plotting.py` and evaluate Numerov's method error (Exercise 5.1-e)[[1]](#1) with `51e.py`.

## References
<a id="1">[1]</a> 
Thijssen, Jos (2007). Computational Physics. 2ª ed. Cambridge University Press. DOI:10.1017/CBO9781139171397.43
