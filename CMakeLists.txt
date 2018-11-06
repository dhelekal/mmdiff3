cmake_minimum_required(VERSION 3.12)
project(mmdiff3)

set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
set(OpenMP_C_FLAGS "-Xclang -fopenmp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

find_package(OpenMP REQUIRED)

if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP_TARGET INTERFACE)
    add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(mmdiff3 SHARED main.cpp mmdiff_3.cpp mmdiff_3.hpp)
add_executable(main main.cpp)

target_link_libraries(mmdiff3 PUBLIC OpenMP::OpenMP_CXX Eigen3::Eigen)