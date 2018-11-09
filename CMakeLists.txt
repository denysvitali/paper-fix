project(paperfix)
cmake_minimum_required(VERSION 3.12)

find_package(OpenCV REQUIRED)

add_executable(paperfix main.cpp)
target_link_libraries(paperfix ${OpenCV_LIBS})
