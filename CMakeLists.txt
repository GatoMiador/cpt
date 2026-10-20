cmake_minimum_required(VERSION 3.18...3.31)

project(cpt LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)

execute_process(
	COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
	OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Python module with the C++ implementation (cpt.h)
nanobind_add_module(_cpt NB_STATIC python/_cpt.cpp)

target_include_directories(_cpt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# M_PI on MSVC
if (MSVC)
	target_compile_definitions(_cpt PRIVATE _USE_MATH_DEFINES)
endif()

install(TARGETS _cpt LIBRARY DESTINATION .)
install(FILES cpt.py DESTINATION .)
