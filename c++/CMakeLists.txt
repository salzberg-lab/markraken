cmake_minimum_required(VERSION 2.8.4)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
project(markraken)

# set(CMAKE_CXX_STANDARD 11)

add_executable(markraken main.cpp HPC.cpp HPC.h include/FastaTools.h include/FastaTools.cpp include/NcbiTaxonomy.h include/NcbiTaxonomy.cpp markerizer.cpp markerizer.h index.cpp index.h)