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
elif [ "$1" == 'set' ]
then
	version=`date +%Y%m%d`
	sed -i CMakeLists.txt \
		-e "s|\(SET( *implicit_ad_version\)* \"[0-9]\{8\}\"|\1 \"$version\"|"
	echo 'bin/version.sh: OK'
	exit 0
elif [ "$1" == 'copy' ]
then
	version=`bin/version.sh get`
	sed -i implicit_ad.omh \
		-e "s|implicit_ad-[0-9]\\{8\\}|implicit_ad-$version|"
	echo 'bin/version.sh: OK'
	exit 0
else
	echo 'usage: bin/version.sh (get|set|copy)'
	exit 1
fi
