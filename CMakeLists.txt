#
# print_variable(variable)
MACRO(print_variable variable)
	MESSAGE(STATUS "${variable} = ${${variable}}" )
ENDMACRO(print_variable)
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(newton_step)
SET(newton_step_version "20171014")
#
# must have pkg-config
FIND_PACKAGE(PkgConfig)
#
# use pkg-config to find include directories
FOREACH(package cppad eigen3)
	pkg_check_modules(${package} REQUIRED ${package})
	print_variable(${package}_INCLUDE_DIRS)
	INCLUDE_DIRECTORIES( ${${package}_INCLUDE_DIRS} )
ENDFOREACH(package)
#
#
ADD_EXECUTABLE(check
	src/check.cpp
	src/control.cpp
	src/utility.cpp
	src/implicit_kedem.cpp
	src/implicit_newton.cpp
)
ADD_EXECUTABLE(time
	src/time.cpp
)
