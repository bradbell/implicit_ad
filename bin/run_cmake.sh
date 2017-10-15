#! /bin/bash -e
# -----------------------------------------------------------------------------
# value used for CMAKE_VERBOSE_MAKEFILE
cmake_verbose_makefile='false'
# value used for CMAKE_BULD_TYPE
cmake_build_type='release'
# extra flags, besides debug and release flags, used during compliation
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
source_dir='..'
cmake \
	-D CMAKE_VERBOSE_MAKEFILE="${cmake_verbose_makefile}" \
	-D CMAKE_BUILD_TYPE="${cmake_build_type}" \
	-D extra_cxx_flags="${extra_cxx_flags}" \
	${source_dir}
