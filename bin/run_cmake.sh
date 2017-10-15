#! /bin/bash -e
# -----------------------------------------------------------------------------
extra_cxx_flags='-Wall -pedantic-errors -std=c++11 -Wshadow'
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/run_cmake.sh' ]
then
	echo 'bin/run_cmake.sh: must be run from its parent directory'
	exit 1
fi
if [ ! -e 'build' ]
then
	mkdir build
fi
cd build
cmake -D extra_cxx_flags="${extra_cxx_flags}" ..
