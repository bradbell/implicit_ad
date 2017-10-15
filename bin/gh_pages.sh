#! /bin/bash -e
# -----------------------------------------------------------------------------
# check usage
if [ "$0" != 'bin/gh_pages.sh' ]
then
	echo 'bin/gh_pages.sh: must be run from its parent directory'
	exit 1
fi
#
# make sure that gh_pages.sh is the only file that has changed
list=`git status -s | sed -e '/ bin\/gh_pages.sh$/d'`
if [ "$list" != '' ]
then
	git status -s | sed -e '/ bin\/gh_pages.sh$/d'
	echo 'gh_pages.sh: git files other than gh_pages.sh have changed'
	exit 1
fi
#
# create an empty build/doc
if [ -e 'build/doc' ]
then
	rm -r build/doc
fi
mkdir build/doc
#
# build copy of current documentation in build/doc
cd build/doc
omhelp ../../implicit_ad.omh -debug -noframe
cd ../..
#
# copy current gh_pages.sh to a safe place
cp bin/gh_pages.sh build/gh_pages.sh
#
echo 'gh_pages.sh: OK'
exit 0
