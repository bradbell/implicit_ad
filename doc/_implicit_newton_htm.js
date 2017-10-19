var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'implicit_ad.htm',
'implicit_newton.htm'
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
var list_down0 = [
'test_circle_implicit_newton.htm'
];
var list_current0 = [
'implicit_newton.htm#Syntax',
'implicit_newton.htm#Purpose',
'implicit_newton.htm#full_step',
'implicit_newton.htm#num_step',
'implicit_newton.htm#aL_fun',
'implicit_newton.htm#F_fun',
'implicit_newton.htm#solve',
'implicit_newton.htm#solve.solve.function',
'implicit_newton.htm#solve.solve.derivative',
'implicit_newton.htm#solve.solve.linear',
'implicit_newton.htm#k',
'implicit_newton.htm#xk',
'implicit_newton.htm#rk',
'implicit_newton.htm#q',
'implicit_newton.htm#w',
'implicit_newton.htm#dw',
'implicit_newton.htm#Example'
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
