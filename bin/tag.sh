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
# check that gh-pages branch is up to date
bin/gh_pages.sh
list=`git status -s`
if [ "$list" != '' ]
then
	echo 'tag.sh: git status -s is not empty (for gh-pages branch)'
	echo 'You must commit or abort changes before creating this tag'
	exit 1
fi
#
# check that chages have been pushed
for branch in gh-pages master
do
	git checkout $branch
	local_hash=`git show-ref $branch | grep 'refs/heads/' | \
		sed -e 's| *refs/heads/.*||'`
	remote_hash=`git show-ref $branch | grep 'refs/remotes/' | \
		sed -e 's| *refs/remotes/.*||'`
	if [ "$local_hash" != "$remote_hash"  ]
	then
		echo "Local and remote version of $branch are different"
		git show-ref $branch
		echo "you must first push or abort your changes."
		exit 1
	fi
done
#
# check if tag already exists
set +e
list=`git tag --list | grep "$version"`
set -e
if [ "$list" != '' ]
then
	echo 'An tags for this version already exists. Delete old tags ?'
	for tag in $list
	do
		echo "	git tag -d $tag"
		echo "	git push --delete origin $tag"
	done
	exit 1
fi
#
# create this tag for master
git tag -a -m "create tag for source code branch" $version master
git push origin $version
#
# create this tag for gh-pages
git tag -a -m "create tag for documentation branch" doc-$version gh-pages
git push origin doc-$version
#
echo 'tag.sh: OK'
exit 0
