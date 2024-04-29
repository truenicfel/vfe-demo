# Variational Feature Extraction Demo

## About

This demo contains an isoline extraction on the earth mantle dataset and a critical
line extraction in the ocean data set as shown in our paper.

## Prerequisites

You need to have a compiled version of [VTK](https://vtk.org/) available. This was
tested with VTK 9.2. 

## Build

1. Create a build directory. You might want to use different ones for Release/Debug.
```
mkdir build
```
2. Run CMake configure and generate. During this you will probably be warned that VTK
was not found. Supply the build folder of your VTK build as a variable to CMake.
```
cd build
cmake ..
```
3. Build it!
```
cmake-- build .
```

## Run

Execute the `EarthMantle` executable or the `Ocean` executable like this while supplying 
earthMantle.vti, ocean.vti and/or oceanSeeds.vtp which can be found in the supplementary material of our paper.
When running you might encounter an error message that certain VTK libraries were not
found. You need to make them available to your executable by adding their path to your
environment.
```
./EarthMantle [PATH_TO_EARTH_MANTLE_VTI]
```
```
./Ocean [PATH_TO_EARTH_OCEAN_VTI] [PATH_TO_OCEAN_SEEDS_VTP]
```

After running you should find an `isoLine.vtp` file and/or a `criticalLines.vtp` file
depending on which of the two executables you executed.