var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'implicit_ad.htm',
'control.htm',
'full_newton.htm'
];
var list_down2 = [
'license.htm',
'run_cmake.sh.htm',
'utility.htm',
'implicit_kedem.htm',
'implicit_newton.htm',
'control.htm',
'time.htm'
];
var list_down1 = [
'vector_matrix.htm',
'objective.htm',
'rec_objective.htm',
'constraint.htm',
'rec_constraint.htm',
'full_newton.htm',
'implicit_solver.htm',
'test_control_reduced_objective.htm'
];
var list_current0 = [
'full_newton.htm#Syntax',
'full_newton.htm#Prototype',
'full_newton.htm#Purpose',
'full_newton.htm#xy_in',
'full_newton.htm#xy_out',
'full_newton.htm#L_fun',
'full_newton.htm#criteria',
'full_newton.htm#max_itr',
'full_newton.htm#L_y',
'full_newton.htm#work',
'full_newton.htm#num_itr',
'full_newton.htm#Example'
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
