# SIMPLE

This is an implementation of the average positronium lifetime image reconstruction method SIMPLE (Statistical IMage reconstruction of Positron Lifetime via time-wEighting).

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [License](#license)
5. [Contact](#contact)

## Overview

This project performs triples pairing and lifetime image reconstruction, whose code locate in ./processing and ./recon respectively.  

## Installation

No installation is required. However, re-make of the mex functions may be necessary depending on your environment. To do this, run the scripts with name 'make*.m' (beginning with make) under folder ./recon/funcs and ./processing/mex and funcstions for list-mode projection and data processing shall be ready.

## Usage

The data processing and image reconstruction are two standalone parts. The later reconstruction part can be executed without running the first part of processing as the list-mode data is supplied.

### Source data

https://ucdavis.box.com/s/hnrpxki5rb4wo1kp85dqu70a642u1r3f

### Data processing



### SIMPLE reconstruction

In ./recon, run_moby.m and run_prismpet.m perform the image reconstruction 

- Download the ./simulation and ./prism-pet folders in the above link to your drive
- Change the 'data_folder' variable to your data paths in run_moby.m and run_prismpet.m respectively
- Hit run

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Acurbbb/SIMPLE/blob/main/LICENSE) file for details.

## Contact

For any questions or suggestions, feel free to contact Bangyan Huang at bybhuang@ucdavis.edu or Jinyi Qi at qi@ucdavis.edu
