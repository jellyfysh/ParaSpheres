cmake_minimum_required(VERSION 3.16)
project(GraphValidateCellECMC)
include_directories(/opt/local/include)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(GraphValidateCellECMC
        cell.cc
        checker.cc
        Checker.h
        GraphValidateCellECMC.cc
        init.cc
        inout.cc
        main.cc
        Param.h
        refresh.cc)
