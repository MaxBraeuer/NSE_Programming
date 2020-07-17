cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# set default compiler
option(USE_CLANG "build application with clang" OFF) # OFF is the default
if(USE_CLANG)
    SET (CMAKE_C_COMPILER   "/usr/bin/clang")
endif(USE_CLANG)

SET (CMAKE_C_FLAGS  "-Wall -Wextra -Os -lm")


project(NAVIER_STOKES_EQUATIONS)

#----------------------------------------------------------------------------
# Use most-recent available language dialect
#
set(CMAKE_C_STANDARD 11)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include) 
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable
#
add_executable(nse NSE.c ${sources} ${headers})

#----------------------------------------------------------------------------
# Copy the UserInput.txt into the build directory
#
set(USER_INPUT UserInput.txt)
configure_file(
    ${PROJECT_SOURCE_DIR}/${USER_INPUT}
    ${PROJECT_BINARY_DIR}/${USER_INPUT}
    COPYONLY)
