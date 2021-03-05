# studious\_potato

This repository contains code developed during my master's thesis.
It does not include Cyclus, Cycamore or Misonenrichment, which are avaible
here (in brackets the version used):

- [Cyclus](https://github.com/cyclus/cyclus) (v1.5.5)
- [Cycamore](https://github.com/cyclus/cycamore) (v1.5.5)
- [Misoenrichment](https://github.com/maxschalz/miso_enrichment) (v1.0)

All Python scripts are written for Python 3.6. Additional Python libraries
are needed such as SciPy, NumPy, Matplotlib, Pandas and Sqlite3.

## data
The `data` folder contains the SERPENT output files featuring the reactor
simulation, as well as digitised data from publications to counter-check 
data.
It does _not_ include the Cyclus output files. However, each folder 
contains the input files, thus the simulations can easily be recreated.
In order for the simulation to work, it must be called from within the 
`input_files` directory.

## analysis
This `folder` contains all the analysis scripts.

## License
The source code, i.e. the code found in `analysis` and `data`, is licensed under the 
[BSD-3 Clause License](https://github.com/maxschalz/studious_potato/blob/main/LICENSE).
The thesis, the defence presentation and the corresponding work, all found in `thesis` and `presentation`, are licensed under the 
[Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
