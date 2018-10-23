cmake_minimum_required(VERSION 3.12)
project(mmdiff3)

set(CMAKE_CXX_STANDARD 17)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(mmdiff3 main.cpp mmdiff_3.cpp mmdiff_3.hpp)

target_link_libraries (mmdiff3 Eigen3::Eigen)