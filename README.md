# Fortran code for the formation of limestone-marl alternations

This is the Fortran code that was used to produce the results in L'Heureux 2018 paper "Diagenetic Self-Organization and Stochastic Resonance in a Model of Limestone-Marl Sequences".

The code has been modified by Lukas van de Wiel (Utrecht University) to match the Fortran 2008 standard. That version is in the `lheureux.f90` file. The `.f90` extension seems to be required by Ifort so it is kept here.

## Building and running

Run `make` to build the code with `gfortran`. The executable should appear in the `./build` directory. Running it by typing `./build/marlpde08` will produce several output files in plain text format. `make` builds the Fortran 2008 version. It gives the same results as the original code, but is slightly faster. 

Run `make clean` to remove all files produced during compilation.

## Citation

If you use this code in your work, please cite both the original paper by Ivan L'Heureux and this repository, following the information in CITATION.cff.

## License

The maintainers of this package are not directly affiliated with the original author, but obtained permission to distribute this code under the Apache 2 license.

>  Copyright 2022 Ivan L'Heureux (University of Ottawa)
>
>   Licensed under the Apache License, Version 2.0 (the "License");
>   you may not use this file except in compliance with the License.
>   You may obtain a copy of the License at
>
>       http://www.apache.org/licenses/LICENSE-2.0
>
>   Unless required by applicable law or agreed to in writing, software
>   distributed under the License is distributed on an "AS IS" BASIS,
>   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
>   See the License for the specific language governing permissions and
>   limitations under the License.

## Code versions

- The original code as sent to as by Ivan L'Heureux is in the protected branch `corrected_archive`. `corrected` refers to corrections made by I. L'Heureux himself upon our request. This is the basis for the code in the `main` branch.

- Johan Hidding expanded the code to save output into a `hdf5` file and allow reading parameters from a config file. However, this seems to cause problems due to the conversion of data types and most likely needs testing and modifications before it can be relied upon. This version is in the `config_file` branch. It includes the `cfg` parser from [pkgpl/cfgio](https://github.com/pkgpl/cfgio) by Wansoo Ha.