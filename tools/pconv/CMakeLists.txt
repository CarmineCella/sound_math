project (pconv)

SET(CMAKE_CXX_FLAGS "-Wall -g2 -O2 -std=c++11")

include_directories("../../include")

add_executable(pconv pconv.cpp)

add_custom_command(TARGET pconv POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pconv> ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/../../LICENSE ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/../../samples/anechoic1.wav ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/../../samples/Concertgebouw-s.wav ${PROJECT_SOURCE_DIR}/../../build/bin/
)
