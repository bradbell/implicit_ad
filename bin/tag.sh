#! /bin/bash -e
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/tag.sh' ]
then
	echo 'bin/tag.sh: must be run from its parent directory'
	exit 1
fi
#
branch=`git branch | sed -n -e '/^\*/p' | sed -e 's|^\* *||'`
if [ "$branch" != 'master' ]
then
	echo 'tag.sh: can only be executed using the master branch'
	exit 1
fi
#
# make sure version is consistent
version=`bin/version.sh get`
bin/version.sh copy $version
#
# make sure there are no uncommited changes
list=`git status -s`
if [ "$list" != '' ]
then
	echo 'tag.sh: git status -s is not empty (for master branch)'
	echo 'You must commit or abort changes before creating this tag'
	exit 1
fi
#
# check if tag already exists
if git tag --list | grep "$version"
then
	echo 'This git referece tag already exists. Delete old version ?'
	echo "	git tag -d $version"
	echo "	git push --delete origin $version"
	exit 1
fi
#
# create this tag
git tag -a -m "create new release" $version
git push origin $version
#
echo 'tag.sh: OK'
exit 0
