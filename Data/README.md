### Data files

This directory contains all the data files produced by the code in this repository and presented in the thesis. 

- `200 GeV`: Full and EPA calculations at $\sqrt{s}=200 \, \mathrm{GeV}$
- `13 TeV`: Full and EPA calculations at $\sqrt{s}=13 \, \mathrm{TeV}$
- `13 TeV (cuts)`: Full calculation at $\sqrt{s}=13 \, \mathrm{TeV}$ with ATLAS cuts

The folders `200 GeV` and `13 TeV` share a similar structure in the files. The full calculation can be found in a file with name `full_me_*.csv` and the EPA calculations can be found in files with name `epa_*_i.csv`, where `i = 1, 2, 3`. The different values of `i` correspond to different photon distributions:

- `i = 1`: Eq. (6.6) with $M=0$, $\mu_p=1$
- `i = 2`: Eq. (6.5) with $\mu_p=1$
- `i = 3`: Eq. (6. 34)

The folder `13 TeV (cuts)` contains four files of the form `full_me_13000_i.csv`, where `i = 1, 2, 3, 4`. Here `i` corresponds to the invariant-mass bin in the ATLAS data.

### CSV file format

All files share a common CSV format, with the headers

[invariant mass bin lower edge], [invariant mass bin upper edge], [cross section in $\mathrm{pb}/\mathrm{GeV}$], [absolute error in $\mathrm{pb}/\mathrm{GeV}$], [$\chi^2$], [number of rebinning iterations]