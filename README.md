# SIMPLE

This is an implementation of the average positronium lifetime image reconstruction method SIMPLE (Statistical IMage reconstruction of Positron Lifetime via time-wEighting).

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [License](#license)
5. [Contact](#contact)

## Overview

This project performs triples pairing and lifetime image reconstruction. Executable scripts locate in ./processing and ./recon respectively.  

## Installation

No installation is required. However, re-make of the mex functions may be necessary depending on your environment. To do this, run the scripts with name 'make*.m' (beginning with make) under folder ./recon/funcs and ./processing/mex and funcstions for list-mode projection and data processing shall be ready.

## Usage

The data processing and image reconstruction are in two standalone parts. The later reconstruction part can be executed without running the first part of processing as the list-mode data is supplied.

### Source data

https://ucdavis.box.com/s/hnrpxki5rb4wo1kp85dqu70a642u1r3f

### Data processing

The processing pipeline can be found in ./processing. 

- Download the ./simulation_raw folder using the above link. The raw singles data is in list-mode with each event being coded in a custom format
- (Optional) Re-make the mex functions for singles data reading and triples pairing
- In singles2triple.m, change the input data folder to your data location and specifiy a output path. Hit Run
- For each of the small raw singles file in the folder 'raw_singles', the previous step generates a .lm file to represent the double in a triple and a .float file for the lifetime measurement. The .lm file has crystal locations (2 transaxial IDs and 2 axial IDs) and TOF info for each double. Each event is represented by 5 int16 numbers. The correponding lifetime measurement is stored in the .float file in the same order
- Change the input path to the location of your .lm and .float files. Run to_final_lm_simple_moment.m. This step generates reconstruction-ready lifetime events in both reconstruction and correction time windows.

### SIMPLE reconstruction

In ./recon, run_moby.m and run_prismpet.m perform the image reconstruction. Each reconstruction script requires three data sources: prompts, delayeds, and total doubles events.

- Download the ./simulation and ./prism-pet folders in the above link to your drive
- (Optional) Re-make the mex functions for list-mode projection
- Change the 'data_folder' variable to your data paths in run_moby.m and run_prismpet.m respectively
- Hit Run

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Acurbbb/SIMPLE/blob/main/LICENSE) file for details.

## Contact

For any questions or suggestions, feel free to contact Bangyan Huang at bybhuang@ucdavis.edu or Jinyi Qi at qi@ucdavis.edu
