cmake_minimum_required(VERSION 3.10)
project(MatrixOperations)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_library(MatrixOperations src/Matrix.cpp)

add_executable(MatrixApp src/main.cpp)
target_link_libraries(MatrixApp MatrixOperations)
