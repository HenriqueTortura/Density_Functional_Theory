# Density Functional Theory

## Build Instructions
Compile main module for Fortran (in my case, f95 calls for GNU Fortran)
```bash
f95 -c dft.f95
```
and the module for Python
```bash
f2py -c -m pydft dft.f95
```
## Running All Problems From Thijssen
To see all problems from Thijssen[[1]](#1), simply run the script
```bash
python3 Solutions.py
```
Each problem (and the section 5.5) can be run after a y/n command.

## Running Specific Cases
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

## References
<a id="1">[1]</a> 
Thijssen, Jos (2007). Computational Physics. 2ª ed. Cambridge University Press. DOI:10.1017/CBO9781139171397.43
