#! /bin/bash -e
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/run_omhelp.sh' ]
then
	echo 'bin/run_omhelp.sh: must be run from its parent directory'
	exit 1
fi
if [ ! -e 'doc' ]
then
	mkdir doc
fi
cd doc
omhelp ../implicit_ad.omh -debug -noframe -xml
omhelp ../implicit_ad.omh -debug -noframe
