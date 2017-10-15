var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'implicit_ad.htm',
'utility.htm',
'solve_lower_cppad.htm'
];
var list_down2 = [
'run_cmake.sh.htm',
'utility.htm',
'implicit_kedem.htm',
'implicit_newton.htm',
'control.htm',
'time.htm'
];
var list_down1 = [
'notation.htm',
'norm_squared.htm',
'join_vector.htm',
'sparse_cppad2eigen.htm',
'solve_lower_cppad.htm',
'jac_constraint.htm'
];
var list_current0 = [
'solve_lower_cppad.htm#Syntax',
'solve_lower_cppad.htm#Prototype',
'solve_lower_cppad.htm#A',
'solve_lower_cppad.htm#x',
'solve_lower_cppad.htm#b',
'solve_lower_cppad.htm#msg',
'solve_lower_cppad.htm#Example'
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
function choose_down2(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down2[index-1];
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
