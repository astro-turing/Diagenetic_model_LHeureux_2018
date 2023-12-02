# L'Heureux Legacy Fortran code
This is the Fortran code that was used to produce the results in L'Heureux 2018 paper "Diagenetic Self-Organization and Stochastic Resonance in a Model of Limestone-Marl Sequences". The original code is `lheureux.f` and it has comments added by Niklas Hohmann.

It has been expanded by Johan Hidding to save output into a `hdf5` file and allow reading parameters from a config file. However, this seems to cause problems due to the conversion of data types.

The code has been modified by Lukas van de Wiel to match the Fortran 2008 standard. That version is in the `lheureux.f90` file. The `.f90` extension seems to be required by Ifort so it is kept here.

## Building and running
Run `make` to build the code with `gfortran`. The executable should appear in the `./build` directory. Running it will produce several output files in plain text format.

`make` builds the Fortran 2008 version without the extensions. Run `make all` for the original version of the code with the config file reader.

## Citation
If you use this code in your work, please cite both the original paper by Ivan L'Heureux and this repository, following the information in CITATION.cff.

## License
The maintainers of this package are not directly affiliated with the original author, but obtained permission to distribute this code under the Apache 2 license.

>  Copyright 2022 Ivan L'Heureux (University of Ottawa), Johan Hidding (Netherlands eScience Center)
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

## 3rd party code
This repository includes the `cfg` parser from [pkgpl/cfgio](https://github.com/pkgpl/cfgio) by Wansoo Ha.
