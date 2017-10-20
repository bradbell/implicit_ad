#! /bin/bash
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
#	cppad
#	config
#	pkgconfig
#	eigen
#	ipopt
# &&
#
# &section Run CMake to Configure Implicit AD&&
#
# &head Syntax&&
# &codei%bin/run_cmake.sh%&&
#
# &head Purpose&&
# This is a
# &href%https://en.wikipedia.org/wiki/Bash_(Unix_shell)%bash%&&
# script that configures &code implicit_ad&&
# by running &href%https://en.wikipedia.org/wiki/CMake%CMake%&&.
# If you do not have bash on your system,
# this script can be considered documentation for the &code cmake&& command.
#
# &subhead Correctness Testing&&
# After this script is run, the following commands will run the
# tests for correctness:
# &codep
#	cd build
#	make check
# &&
#
# &subhead Speed Testing&&
# After this script is run, the following commands will run the
# speed tests:
# &codep
#	cd build
#	make time
#	./time
# &&
#
# &head Exit on Error&&
# The following command instructs bash to terminate, with a non-zero
# exit status, when an error occurs:
# &srccode%sh%
set -e
# %&&
#
# &head Check Working Directory&&
# Check that the current working directory is the top level source directory.
# &srccode%sh%
if [ "$0" != 'bin/run_cmake.sh' ]
then
	echo 'bin/run_cmake.sh: must be run from its parent directory'
	exit 1
fi
# %&&
#
# &head cppad_pkg_config_path&&
# This is the directory where the file &code cppad.pc&& file is located
# &srccode%sh%
cppad_pkg_config_path="$HOME/prefix/cppad/share/pkgconfig"
# %&&
#
# &head eigen_pkg_config_path&&
# This is the directory where the file &code eigen3.pc&& file is located.
# It is often desirable to install eigen its own special prefix
# so warnings can be suppressed for its include files without suppressing
# warnings for other include files.
# &srccode%sh%
eigen_pkg_config_path="$HOME/prefix/eigen/share/pkgconfig"
# %&&
#
# &head ipopt_pkg_config_path&&
# The install of &code cppad&& may be done in a way that uses the
# &code ipopt&& libraries. If this is the case, you must change the following
# setting to the directory where the file &code ipopt.pc&& is located
# (otherwise, the &code cmake&& command will report the error message
# &code 'ipopt' required by 'cppad'&&).
# &srccode%sh%
ipopt_pkg_config_path=""
# %&&
#
# &head PKG_CONFIG_PATH&&
# Set the directories that are searched by
# &href%https://en.wikipedia.org/wiki/Pkg-config%pkg-config%&&.
# &srccode%sh%
PKG_CONFIG_PATH="${cppad_pkg_config_path}:${eigen_pkg_config_path}"
if [ "$ipopt_pkg_config_path" != "" ]
then
	PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${ipopt_pkg_config_path}"
fi
# %&&
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
cmake_build_type='debug'
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
# &head Exit &&
# This script runs with
# &srccode%sh%
echo 'run_cmake.sh: OK'
exit 0
# %&&
#
# &end
