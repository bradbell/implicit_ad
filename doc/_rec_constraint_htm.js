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
'rec_constraint.htm'
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
'rec_constraint.htm#Syntax',
'rec_constraint.htm#Prototype',
'rec_constraint.htm#L_fun',
'rec_constraint.htm#L_fun.J',
'rec_constraint.htm#L_fun.delta_t',
'rec_constraint.htm#p',
'rec_constraint.htm#p.x',
'rec_constraint.htm#p.y',
'rec_constraint.htm#p.L',
'rec_constraint.htm#Example',
'rec_constraint.htm#Example.ADFun&lt;double&gt;',
'rec_constraint.htm#Example.ADFun&lt; AD&lt;double&gt; &gt;'
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
