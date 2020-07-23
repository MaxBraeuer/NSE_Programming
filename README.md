# Numerical Solution of the Navier-Stokes Equations

*Einleitende Beschreibung*

## Building and Running

Uses CMake for out-of-place builds with a build environment of your choice. This is an example work-flow and configuration to build the NSE source:

1. From the source directory create a build directory, enter it and additionally create a directory called `output` where the program is going to save its raw data:

    * `mkdir build && cd build && mkdir output`

2. Run CMake: 

    * `cmake -G "<generator>" [options] ..` 

    Some common build system generators are:

    * `Ninja` --- for generating Ninja build files.
    * `Unix Makefiles` --- for generating make-compatible parallel makefiles.
    * `Visual Studio` --- for generating Visual Studio projects and solutions.
    * `Xcode` --- for generating Xcode projects.

    To change the default C-compiler or to set additional compiler flags (optionally) edit the `CMakeLists.txt` file in the source directory or use command line options.
    For possible options consult the [CMake documentation](https://cmake.org/documentation/) or execute `cmake --help`.

3. Build the application by invoking your build system specified above directly, e.g.

    * `ninja` --- Ninja will run a parallel build by default.
    * `make` --- for a parallel build with `NN` parallel jobs execute `make -j NN`

4. To run the application simply execute `./nse` from within your build directory.

## Settings and Parameters

User configurable variables can be set post-build by editing the file `UserInput.txt` *located within the build-directory.*

## Results

...
