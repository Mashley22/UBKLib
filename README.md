# UBKLib
A library for efficiently tracing the Earth's magnetic field lines and 
calculating the second adiabatic invariant along the field line. The final 
intention is to be able to use the method described in [Whipple 1978](https://doi.org/10.1029/JA083iA09p04318)
to be able to find close paths of particles in the Earth's magnetosphere.

# Description
The core tracer is written using c++ 23, making use of c++ modules. 
The tracer uses templates to do much configuration, guaranteeing compile time
optimisations. The method is fairly easy and safely parrallelised (TODO) and is 
quite efficient. 

# Field Model Support
The tracer currently has in build support for the [IGRF13](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field)
and [Tsyganenko89](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/t89/)
field models.

# Usage
## Building
Due to the c++ 23 requirements, and the use of c++ 23 modules, the latest possible
compiler versions are recommended. Cmake is used for building, and as such make is
not currently a valid generator (due to the c++ modules), use ninja instead. 
The build directory provides an easy run bash script to build unit tests and samples.
It will also attempt to generate compile_commands.json if clang is available, also
prefferring to use clang if available.

The use of O3 and lto should be stable, no-math-errno is strongly recommended, 
offering around a 10x speedup on the igrf13 model specifically. No testing has 
been done using fast-math.

# Dependencies
All dependencies are added via git submodules. Currently:
 - [Catch2](https://github.com/catchorg/Catch2) : unit testing
 - [cxform2](https://github.com/Mashley22/cxform2) : for some coordinate transforms
 - [nessan/xoshiro](https://github.com/nessan/xoshiro) : for efficient random number generation
