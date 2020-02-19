cmake_minimum_required(VERSION 3.10.2)

project(S111-Streamlines)

find_package(HDF5 REQUIRED COMPONENTS CXX)

include_directories(${HDF5_INCLUDE_DIRS})

include_directories(include)

# MSVC requires following definition before including cmath to define M_PI and such
add_definitions(-D_USE_MATH_DEFINES)

add_executable(s111_to_streamlines src/s111_to_streamlines.cpp)

target_link_libraries(s111_to_streamlines ${HDF5_CXX_LIBRARIES})
