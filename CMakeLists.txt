#
# preamble
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(implicit_ad)
SET(implicit_ad_version "20171014")
#
# find pkg-config
FIND_PACKAGE(PkgConfig)
#
# use pkg-config to find include directories
FOREACH(package cppad eigen3)
	pkg_check_modules(${package} REQUIRED ${package})
	MESSAGE(STATUS "${package}_INCLUDE_DIRS = ${${package}_INCLUDE_DIRS}")
	INCLUDE_DIRECTORIES( ${${package}_INCLUDE_DIRS} )
ENDFOREACH(package)
#
# check
ADD_EXECUTABLE(check
	src/check.cpp
	src/control.cpp
	src/utility.cpp
	src/implicit_kedem.cpp
	src/implicit_newton.cpp
)
# time
ADD_EXECUTABLE(time
	src/time.cpp
)
