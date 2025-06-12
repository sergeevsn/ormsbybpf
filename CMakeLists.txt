# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.14)

# Project name
project(MYPROJECT)

# Set C++ standard (e.g., C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable optimization -O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find FFTW3 library
find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3)

if(NOT FFTW3_LIBRARY)
    message(FATAL_ERROR "FFTW3 library not found!")
endif()

# Build executable from source files
add_executable(test src/test.cpp src/ormsbybpf.cpp)

# Link FFTW3 and math library (required by FFTW)
target_link_libraries(test PRIVATE ${FFTW3_LIBRARY} m)