cmake_minimum_required(VERSION 3.10) # specify the minimum version of CMake that your project requires
project(integer-complexity) # name your project

set(CMAKE_CXX_STANDARD 11) # specify the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(msieve STATIC IMPORTED)
set_target_properties(msieve PROPERTIES
        IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/msieve/libmsieve.a
        INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/msieve/include
)

include_directories(msieve/include) # specify the directories where the compiler should look for include files

add_executable(exp exp.cpp experiment.cpp exp_container.cpp fretriever.cpp helpers.cpp factorize.cpp) # specify the source files for your executable

target_link_libraries(exp gmpxx gmp msieve z ecm pthread png dl) # specify the libraries to link against