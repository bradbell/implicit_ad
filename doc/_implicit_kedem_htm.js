var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'implicit_ad.htm',
'implicit_kedem.htm'
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
'test_circle_implicit_kedem.htm'
];
var list_current0 = [
'implicit_kedem.htm#Syntax',
'implicit_kedem.htm#Purpose',
'implicit_kedem.htm#L_fun',
'implicit_kedem.htm#F_fun',
'implicit_kedem.htm#solve',
'implicit_kedem.htm#solve.solve.function',
'implicit_kedem.htm#solve.solve.derivative',
'implicit_kedem.htm#solve.solve.linear',
'implicit_kedem.htm#k',
'implicit_kedem.htm#xk',
'implicit_kedem.htm#store',
'implicit_kedem.htm#rk',
'implicit_kedem.htm#Method',
'implicit_kedem.htm#Method.Zero Order',
'implicit_kedem.htm#Method.Higher Orders',
'implicit_kedem.htm#Example'
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
