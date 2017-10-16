#! /bin/bash -e
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/version.sh' ]
then
	echo 'bin/version.sh: must be run from its parent directory'
	exit 1
fi
if [ "$1" == 'get' ]
then
	sed -n -e '/SET( *implicit_ad_version/p' CMakeLists.txt |
		sed -e 's/SET( *implicit_ad_version *"//' -e 's/" *)//'
	exit 0
fi
if [ "$1" == 'set' ]
then
	version=`date +%Y%m%d`
	sed -i CMakeLists.txt \
		-e "s|\(SET( *implicit_ad_version\)* \"[0-9]\{8\}\"|\1 \"$version\"|"
fi
if [ "$1" == 'copy' ]
then
	version=`bin/version.sh get`
	sed -i implicit_ad.omh \
		-e "s|implicit_ad-[0-9]\\{8\\}|implicit_ad-$version|"
fi
