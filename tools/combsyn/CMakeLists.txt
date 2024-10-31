project (combsyn)

include_directories("../../include")

add_executable(combsyn combsyn.cpp)

add_custom_command(TARGET combsyn POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:combsyn> ${PROJECT_SOURCE_DIR}/../../build/bin/
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/../../LICENSE ${PROJECT_SOURCE_DIR}/../../build/bin/
)
