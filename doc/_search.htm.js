// ------------------------------------------------------------ 
// Copyright (C) Bradley M. Bell 1998-2015, All rights reserved 
// ------------------------------------------------------------ 
Keyword = 
[
'implicit_ad  AD Methods That Differentiate Implicit Functions  ',' git repository current version ',
'license  Implicit AD License  ',' ',
'run_cmake.sh  Run CMake to Configure Implicit AD  ',' syntax purpose correctness testing speed exit error check working directory cppad_pkg_config_path eigen_pkg_config_path ipopt_pkg_config_path cmake_verbose_makefile cmake_build_type extra_cxx_flags usage command ',
'utility  Utilities Used by All Methods  ',' ',
'notation  Notation  ',' dynamic sparsematrix vector cppad_sparse ',
'norm_squared  Norm Squared of a Vector  ',' ',
'join_vector  Join Two Vectors  ',' ',
'sparse_cppad2eigen  Convert A CppAD Sparse Matrix to an Eigen Sparse Matrix  ',' syntax prototype sparse_eigen example ',
'solve_lower_cppad  Solve a CppAD Sparse Lower Triangular System  ',' syntax prototype b msg example ',
'jac_constraint  Compute Jacobian of Implicit Function Constraints  ',' syntax prototype purpose l_y l_fun xy work example simple control problem ',
'implicit_kedem  Kedem Method for Derivatives of Implicit Functions  ',' syntax purpose l_fun f_fun solve solve.function solve.derivative solve.linear xk store rk zero order higher orders example ',
'test_circle_implicit_kedem  Example / Test of Implicit Wagner Class  ',' ',
'implicit_newton  Newton Step Method for Derivatives of Implicit Functions  ',' syntax purpose full_step num_step al_fun f_fun solve solve.function solve.derivative solve.linear k xk rk q dw example ',
'test_circle_implicit_newton  Example / Test of Implicit Newton Class  ',' ',
'control  The Control Test Problem  ',' ',
'vector_matrix  Conversions Between Control Vectors and Matrices  ',' xy_vec2mat ',
'objective  Computes the Control Objective Function  ',' syntax prototype definition example ',
'rec_objective  Record the Control Objective  ',' syntax prototype f_fun delta_t q ',
'constraint  Computes the Control Constraint Function  ',' syntax prototype definition example ',
'rec_constraint  Record the Control Constraint as a CppAD Function Object  ',' syntax prototype l_fun delta_t example adfun<double> ad<double> ',
'full_newton  Execute Full Newton Steps For Control Constraint  ',' syntax prototype purpose xy_in xy_out l_fun criteria max_itr l_y work num_itr example ',
'implicit_solver  Control Problem Solver for Implicit Kedem or Newton Object  ',' syntax purpose l_fun al_fun criteria solve.function solve.derivative solve.linear example ',
'test_control_reduced_objective  Example / Test of Control Problem Reduced Objective  ',' purpose source ',
'time  Timing Comparison of Methods  ',' ',
'set_T_p_and_q  Set T, p, and q  ',' ',
'repeat_kedem_gradient  Repeated Computation of Control Problem Gradient Using Kedem Method  ',' syntax j size ',
'repeat_newton_gradient  Repeated Computation of Control Problem Gradient Using Newton Method  ',' syntax j reverse size ',
'repeat_kedem_hessian  Repeated Computation of Control Problem Hessian Using Kedem Method  ',' syntax j size ',
'repeat_newton_hessian  Repeated Computation of Control Problem Hessian Using Newton Method  ',' syntax j full_step size '
]

var MaxList = 100;
var Nstring = -1;
var Nkeyword = Keyword.length / 2;
Initialize();

function Initialize()
{
	UpdateList();
	document.search.keywords.focus();
}
function UpdateList(event)
{
	key = 0;
	if( window.event )
		key = window.event.keyCode;
	else if( event )
		key = event.which;
	if( key == 13 )
	{	Goto();
		return;
	}
	var string  = document.search.keywords.value;
	if( Nstring == string.length )
		return;
	Nstring     = string.length;

	var word    = string.match(/\S+/g);
	var nword   = 0;
	if(word != null )
		nword   = word.length;

	var pattern = new Array(nword);
	for(var j = 0; j < nword; j++)
		pattern[j] = new RegExp(word[j], 'i');

	var nlist = 0;
	var list  = '';
	for(i = 0; (i < Nkeyword) && (nlist < MaxList) ; i++)
	{
		var match = true;
		for(j = 0; j < nword; j++)
		{	var flag = pattern[j].test(Keyword[2*i]);
			flag     = flag || pattern[j].test(Keyword[2*i+1]);
			match    = match && flag;
		}

		if( match )
		{
			line  = Keyword[2*i].split(/\s+/);
			line  = line.join(' ');
			list  = list + line + '\n';
			nlist = nlist + 1;
		}
	}
	document.search.list.value    = list;
}
function Choose(textarea)
{	var start_select = textarea.value.substring(0, textarea.selectionStart);
	var start_pos    = Math.max(0, start_select.lastIndexOf('\n') );
	var length       = textarea.value.length;
	var end_select   = 
		textarea.value.substring(textarea.selectionEnd, length);
	var end_pos      = end_select.indexOf('\n');
	if( end_pos >= 0 ) 
	{	end_pos = textarea.selectionEnd + end_pos;
	} else 
	{	end_pos = length;
	}
	// highlight the selected line
	textarea.selectionStart = start_pos;
	textarea.selectionEnd   = end_pos;
	// get the choice from the beginning of the line
	var line = textarea.value.substring(start_pos, length);
	var end_choice = line.indexOf(' ');
	if( end_choice >= 0 )
	{	var choice = line.substring(0, end_choice);
		document.search.choice.value = choice.toLowerCase();
	}
	
	return true;
}
function Goto()
{
parent.location = document.search.choice.value + '.htm';
}
