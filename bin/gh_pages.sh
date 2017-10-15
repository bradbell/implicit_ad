#! /bin/bash -e
# Script that updates the gh-pages branch doc directory
#
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# check usage
if [ "$0" != 'bin/gh_pages.sh' ]
then
	echo 'bin/gh_pages.sh: must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# make sure that gh_pages.sh is the only file that has changed
list=`git status -s | sed -e '/ bin\/gh_pages.sh$/d'`
if [ "$list" != '' ]
then
	git status -s | sed -e '/ bin\/gh_pages.sh$/d'
	echo 'gh_pages.sh: git files other than gh_pages.sh have changed'
	exit 1
fi
#
branch=`git branch | sed -n -e '/^\*/p' | sed -e 's|^\* *||'`
if [ "$branch" != 'master' ]
then
	echo 'gh_pages.sh: can only be executed using the master branch'
	exit 1
fi
hash=`git rev-parse HEAD`
#
# create an empty build/doc
if [ -e 'build/doc' ]
then
	rm -r build/doc
fi
echo_eval mkdir -p build/doc
#
# build copy of current documentation in build/doc
cd build/doc
omhelp ../../implicit_ad.omh -debug -noframe
cd ../..
#
# copy current gh_pages.sh to a safe place and then wipe out changes
cp bin/gh_pages.sh build/gh_pages.sh
git reset --hard
#
# checkout gh-pages branch
git checkout gh-pages
#
# determine which files to remove
list=`ls -a doc`
for file in $list
do
	if [ ! -e "build/doc/$file" ]
	then
		echo_eval git rm doc/$file
	fi
done
#
# copy new version of files
list=`ls -a build/doc | sed -e '/^\.*$/d'`
for file in $list
do
	if [ ! -e "doc/$file" ]
	then
		echo "add doc/$file"
	fi
	cp build/doc/$file doc/$file
	git add doc/$file
done
#
#
cat << EOF
The following command can be used to see the changes to gh-pages branch:
	git diff

The following command can be used to commit the chnages to gh-pages branch:
	git commit -m 'update gh-pages to master $hash'

The following commands can be used to return to master and restore gh_pages.sh:
	git checkout master
	cp build/gh_pages.sh bin/gh_pages.sh
EOF
#
echo 'gh_pages.sh: OK'
exit 0
