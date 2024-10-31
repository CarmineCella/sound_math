cmake_minimum_required (VERSION 2.6)

project (pconv)


include_directories("../../include/dsp")
include_directories("../../include/inout")
include_directories("../../include/math")
include_directories("../../include/metatools")
include_directories("../../include/dsp")
include_directories("../../include/inout")

add_executable(pconv pconv.cpp)
