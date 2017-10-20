var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'implicit_ad.htm',
'run_cmake.sh.htm'
];
var list_down1 = [
'license.htm',
'run_cmake.sh.htm',
'utility.htm',
'implicit_kedem.htm',
'implicit_newton.htm',
'control.htm',
'time.htm'
];
var list_current0 = [
'run_cmake.sh.htm#Syntax',
'run_cmake.sh.htm#Purpose',
'run_cmake.sh.htm#Purpose.Correctness Testing',
'run_cmake.sh.htm#Purpose.Speed Testing',
'run_cmake.sh.htm#Exit on Error',
'run_cmake.sh.htm#Check Working Directory',
'run_cmake.sh.htm#cppad_pkg_config_path',
'run_cmake.sh.htm#eigen_pkg_config_path',
'run_cmake.sh.htm#ipopt_pkg_config_path',
'run_cmake.sh.htm#PKG_CONFIG_PATH',
'run_cmake.sh.htm#cmake_verbose_makefile',
'run_cmake.sh.htm#cmake_build_type',
'run_cmake.sh.htm#extra_cxx_flags',
'run_cmake.sh.htm#Check Usage',
'run_cmake.sh.htm#build Directory',
'run_cmake.sh.htm#CMake Command',
'run_cmake.sh.htm#Exit'
];
function choose_across0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_across0[index-1];
}
function choose_up0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_up0[index-1];
}
function choose_down1(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down1[index-1];
}
function choose_down0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down0[index-1];
}
function choose_current0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_current0[index-1];
}
