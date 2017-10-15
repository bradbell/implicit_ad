#
# preamble
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(implicit_ad)
SET(implicit_ad_version "20171014")
#
# print_variable
MACRO(print_variable variable_name)
	MESSAGE(STATUS "${variable_name} = ${${variable_name}}")
ENDMACRO(print_variable)
#
# cmake variables that affect build
print_variable(CMAKE_VERBOSE_MAKEFILE)
print_variable(CMAKE_BUILD_TYPE)
#
SET(extra_cxx_flags "" CACHE STRING
	"Additional C++ compile flags, besides debug and release flags"
)
print_variable(extra_cxx_flags)
#
# find pkg-config
FIND_PACKAGE(PkgConfig)
#
# use pkg-config to find cppad include directories
pkg_check_modules(cppad REQUIRED cppad)
print_variable(cppad_INCLUDE_DIRS)
INCLUDE_DIRECTORIES( ${cppad_INCLUDE_DIRS} )
#
# use pkg-config to find eigen include directories
pkg_check_modules(eigen3 REQUIRED eigen3)
print_variable(eigen3_INCLUDE_DIRS)
INCLUDE_DIRECTORIES( SYSTEM ${eigen3_INCLUDE_DIRS} )
#
# check
ADD_EXECUTABLE(check
	src/check.cpp
	src/control.cpp
	src/utility.cpp
	src/implicit_kedem.cpp
	src/implicit_newton.cpp
)
SET_TARGET_PROPERTIES(check PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}")
#
# time
ADD_EXECUTABLE(time
	src/time.cpp
)
SET_TARGET_PROPERTIES(time PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}")
