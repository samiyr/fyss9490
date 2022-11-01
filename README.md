# Equivalent photon approximation in one- and two-photon exchange processes

This repository contains the relevant code and data files for my Master's thesis on the equivalent photon approximation.

### Code files

The main code files are `FullME.cc` and `EPA.cc`. These calculate the (approximately) differential cross sections $\mathrm{d}\sigma_{\mu^+\mu^-}/\mathrm{d}W_{\mu^+\mu^-}$ for muon pair production in proton-proton collisions using the full matrix element (Eq. (6.31)) and the equivalent photon approximation (Eq. (6.40)), respectively. 

Alongside these main files are a few utility files, which contain code common to both calculations.

The code computes the differential cross section in units $\mathrm{pb}/\mathrm{GeV}$ **without** dividing by bin width. Division by bin width must be done afterwards.

### Dependencies

The integrals are calculated using GSL's VEGAS algorithm, so GSL must be installed. The code has been tested on GSL 2.7.1. 

In addition, `OpenMP` is used for parallelization.

### Building

The files `FullME.cc` and `EPA.cc` must be built separately, but the same build command should work for both. 

On a UNIX system, the following command can be used for building either `FullME.cc` or `EPA.cc`:

```
g++ -std=c++17 -Wall -lgsl -lgslcblas -I/path/to/gsl/include -L/path/to/gsl/lib -fopenmp -lgomp -O3 -ffast-math -DNDEBUG <file>.cc -o <output>
```

where `/path/to/gsl` leads to the GSL installation direction. `<file>.cc` and `<output>` should be replaced by appropriate names. 

On macOS, the flag `-Xpreprocessor` might be necessary, and the flag `-lgomp` should be replaced by `-lomp`.

Removing the flag `-DNDEBUG` enables assertions.