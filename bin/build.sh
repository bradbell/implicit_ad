#! /bin/bash -e
xml='no'     # generate html or xml documentation files
clang='no'   # use clang++ or g++ for compiling
debug='yes'  # build c++ using debug mode
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
list='
	adolc
	check
	inverse_wagner
	example
	full
	half
	time
	wagner
'
match='no'
for name in $list
do
	if [ "$name" == "$1" ]
	then
		match="$name"
	fi
done
if [ "$match" == 'no' ]
then
	echo 'usage: ./build.sh option'
	echo 'where option is one of following:'
	for name in $list
	do
		echo "	$name"
	done
	exit 1
fi
name="$match"
# -----------------------------------------------------------------------------
if [ ! -e build ]
then
	mkdir build
fi
# -----------------------------------------------------------------------------
# C++ cases
if [ "$debug" == 'yes' ]
then
	cxx_flags="-g -Wall -pedantic-errors -std=c++11 -Wshadow"
else
	cxx_flags="-O2 -DNDEBUG -Wall -pedantic-errors -std=c++11 -Wshadow"
fi
if [ "$xml" == 'yes' ]
then
	omhelp_flags='-noframe -xml'
else
	omhelp_flags='-noframe'
fi
if [ "$clang" == 'yes' ]
then
	cxx_compiler='clang++'
else
	cxx_compiler='g++'
fi
prefix="$HOME/prefix"
# -----------------------------------------------------------------------------
# C++ cases
#
if [ "$name" == adolc ]
then
	echo_eval cd build
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include -I$prefix/adolc/include  \
		-L$prefix/adolc/lib64 -l adolc \
		../src/adolc.cpp -o adolc
	echo_eval ./adolc
fi
if [ "$name" == example ]
then
	echo_eval cd build
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include ../src/$name.cpp -o $name
	echo_eval ./$name
fi
if [ "$name" == wagner ]
then
	echo_eval cd build
	if which omhelp >& /dev/null
	then
		omhelp ../src/$name.omh $omhelp_flags
	else
		echo 'Cannot create build/$name.htm because omhelp not in path'
	fi
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include ../src/$name.cpp -o $name
	echo_eval ./$name
fi
if [ "$name" == inverse_wagner ]
then
	echo_eval cd build
	if which omhelp >& /dev/null
	then
		omhelp ../src/$name.hpp $omhelp_flags
	else
		echo 'Cannot create build/$name.htm because omhelp not in path'
	fi
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include \
		-isystem $prefix/eigen/include \
		../src/$name.cpp -o $name
	echo_eval ./$name
fi
if [ "$name" == time ]
then
	echo_eval cd build
	if which omhelp >& /dev/null
	then
		omhelp ../src/$name.omh $omhelp_flags
	else
		echo 'Cannot create build/$name.htm because omhelp not in path'
	fi
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include \
		-isystem $prefix/eigen/include \
		../src/$name.cpp -o $name
	echo_eval ./$name
fi
if [ "$name" == check ]
then
	source_list='
		check.cpp
		control.cpp
		utility.cpp
		implicit_kedem.cpp
		implicit_newton.cpp
	'
	echo_eval cd build
	if which omhelp >& /dev/null
	then
		omhelp ../src/$name.cpp $omhelp_flags
	else
		echo 'Cannot create build/$name.htm because omhelp not in path'
	fi
	src_list=''
	for source in $source_list
	do
		src_list="../src/$source $src_list"
	done
	echo_eval $cxx_compiler $cxx_flags \
		-I$prefix/cppad/include \
		-isystem $prefix/eigen/include \
		$src_list -o $name
	echo_eval ./$name
fi
echo 'build.sh: OK'
exit 0
