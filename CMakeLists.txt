cmake_minimum_required(VERSION 3.2)

project(photoelasticity)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g -Wall")
add_executable(photoelasticity Main.cpp)
set_property(TARGET photoelasticity PROPERTY CXX_STANDARD 17)
target_link_libraries(photoelasticity PUBLIC -Wall -lSDL2 -lpthread)




