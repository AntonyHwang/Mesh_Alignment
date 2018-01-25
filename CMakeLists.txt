cmake_minimum_required(VERSION 3.9)
project(ICP_Alignment)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

#add_executable(ICP_Alignment main.cpp)

add_executable(ICP_Alignment ${SOURCES})