# ########################################################################################
# Use the code in this section for standalone projects
# ########################################################################################
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(MDCsimulation)
find_package(Garfield REQUIRED)
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
# --- Build executables ------------------------------------------------------------
add_executable(MDC MDC.cpp)
target_link_libraries(MDC Garfield::Garfield ${ROOT_LIBRARIES})
add_executable(generate generate.C)
target_link_libraries(generate Garfield::Garfield ${ROOT_LIBRARIES})
# --- copy gas file ------------------------------------------------------------
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Gasfile)
file(GLOB gasfile RELATIVE ${PROJECT_SOURCE_DIR} *.gas)
file(COPY ${gasfile} DESTINATION ${PROJECT_BINARY_DIR}/Gasfile)





