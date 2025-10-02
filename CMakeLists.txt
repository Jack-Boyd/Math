cmake_minimum_required(VERSION 3.10)
project(Math)

set(CMAKE_CXX_STANDARD 23)

include_directories(include)

add_library(Math 
  src/Matrix.cpp
  src/Polynomial.cpp
)

add_executable(MathLibrary src/main.cpp)
target_link_libraries(MathLibrary Math)
