#! /bin/bash -e
# $OMhelpKeyCharacter=&
# &begin run_cmake.sh&& &newlinech #&&
# &spell
#	cmake
#	configures
#	makefile
#	cxx
#	std
#	Wshadow
#	fi
#	mkdir
#	dir
# &&
#
# &section Run CMake to Configure Implicit AD&&
#
# &head Syntax&&
# &codei%bin/run_cmake.sh%&&
#
# &head Purpose&&
# This is a bash shell script that configures &code implicit_ad&&
# by running CMake. If you do not have bash on your system,
# this script can be considered documentation for the &code cmake&& command.
#
# &head cmake_verbose_makefile&&
# If this is &code true&&, make will display the compile and link commands
# that it uses to build the tests.
# &srccode%sh%
cmake_verbose_makefile='false'
# %&&
#
# &head cmake_build_type&&
# You should use &code debug&& for correctness testing and
# &code release&& for speed testing.
# &srccode%sh%
cmake_build_type='release'
# %&&
#
# &head extra_cxx_flags&&
# Extra C++ flags, besides the debug and release flags,
# used during compile commands.
# &srccode%sh%
extra_cxx_flags='-Wall -pedantic-errors -std=c++11 -Wshadow'
# %&&
#
# &head Check Usage&&
# Make sure that &code run_cmake.sh&& is run from its parent directory.
# &srccode%sh%
if [ "$0" != 'bin/run_cmake.sh' ]
then
	echo 'bin/run_cmake.sh: must be run from its parent directory'
	exit 1
fi
# %&&
#
# &head build Directory&&
# If necessary, create the build directory,
# then make it the current working directory.
# &srccode%sh%
if [ ! -e 'build' ]
then
	mkdir build
fi
cd build
# %&&
#
# &head CMake Command&&
# &srccode%sh%
source_dir='..'
cmake \
	-D CMAKE_VERBOSE_MAKEFILE="${cmake_verbose_makefile}" \
	-D CMAKE_BUILD_TYPE="${cmake_build_type}" \
	-D extra_cxx_flags="${extra_cxx_flags}" \
	${source_dir}
# %&&
#
# &end
