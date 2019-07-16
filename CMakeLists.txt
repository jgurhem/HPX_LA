cmake_minimum_required(VERSION 3.0)
project(HPX_LA CXX)

find_package(HPX REQUIRED)
link_directories(${HPX_LIBRARY_DIR})

add_executable(lu_tiled lu_tiled.cpp)
hpx_setup_target(lu_tiled COMPONENT_DEPENDENCIES iostreams)

add_executable(lu_tiled_dist lu_tiled_dist.cpp)
hpx_setup_target(lu_tiled_dist COMPONENT_DEPENDENCIES iostreams)
