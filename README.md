# Garfield++ simulation for MDC

## What is Garfield++

Garfield++ is a toolkit for the detailed simulation of particle detectors based on ionisation measurement in gases or semiconductors. The main area of application is currently in micropattern gaseous detectors.

## What is MDC

The MDC (main drift chamber) is the main part of the tracking system in **STCF**, providing important functions including the following:

- *Reconstructing tracks of charged particles*
- *Measuring energy loss for PID*
- *Providing trigger decision*

## Required dependencies

- [Garfield++](https://garfieldpp.web.cern.ch/garfieldpp/ "Garfield++ official website")
- [ROOT 6](https://root.cern/ "ROOT official website")
- GSL (GNU Scientific Library)
- C++14 or higher
- CMAKE 3.9 or higher

## Function description

### Generate gas file

**generate.C** can generate gas files of different gas mixtures (available gases can be referenced in [Garfield guide](https://garfieldpp.web.cern.ch/garfieldpp/documentation/UserGuide.pdf)) and with many parameters (temperature, pressure, EM field etc.).

### main function

**MDC.cpp** can simulate fundamental element of MDC, (i.e., a drift cell), including geometry definition,primary ionisation, electron/ion transportion, induced signal collection.

## How to run

1. Download or clone source code in Linux system. (suppose absolute path is */home/username/MDCsimulation*)
2. Switch working directory. (`cd /home/username/MDCsimulation`)
3. Make build directory and enter into. (`mkdir build` && `cd ./build`)
4. Cmake and generate executable. (`cmake ..` && `make`)
5. Run the executable (`./MDC gasfile fileid`), you can specify no more than 2 arguments: *gasfile* is the gasfile path, and *fileid* is the suffix of output file.




