project (sounder)

SET(CMAKE_CXX_FLAGS "-Wall -g2 -O2 -std=c++11")

include_directories("../../include")

add_executable(sounder sounder.cpp)
