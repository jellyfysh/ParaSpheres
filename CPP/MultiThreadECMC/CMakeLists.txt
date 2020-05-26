cmake_minimum_required(VERSION 3.13)  # CMake version check
project(simple_example)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 20)            # Enable c++14 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES MultiThreadECMC.cc inout.cc main.cc simul.cc singlethread.cc)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(MultiThreadECMC ${SOURCE_FILES})
