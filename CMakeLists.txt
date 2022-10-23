cmake_minimum_required(VERSION 3.20)
project(labs)

set(CMAKE_CXX_STANDARD 20)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(labs main.cpp)
target_link_libraries (labs Eigen3::Eigen)