/*
$begin check$$
$aindex section head subhead$$
$latex \newcommand{\B}[1]{{\bf #1}}$$
$latex \newcommand{\R}[1]{{\rm #1}}$$

$section Check Implicit Function Derivative Computations$$

$contents%src/utility.hpp
	%implicit_kedem.hpp
	%implicit_newton.hpp
	%control.hpp
	%control.cpp
%$$

$end
*/
// ---------------------------------------------------------------------------
// run tests
// ---------------------------------------------------------------------------
// test runner
# include <cppad/utility/thread_alloc.hpp>
# include <cppad/utility/test_boolofvoid.hpp>

// utility routines
extern bool test_sparse_cppad2eigen(void);
extern bool test_solve_lower_cppad(void);
extern bool test_jac_constraint(void);
extern bool test_control_jac_constraint(void);
//
extern bool test_control_objective(void);
extern bool test_control_constraint(void);
extern bool test_control_rec_constraint(void);
extern bool test_control_full_newton(void);
extern bool test_control_implicit_solver(void);
extern bool test_control_ad_rec_constraint(void);
extern bool test_control_reduced_objective(void);
//
extern bool test_cricle_implicit_kedem(void);
extern bool test_cricle_implicit_newton(void);
//
int main(void)
{
	// create test runner
	CppAD::test_boolofvoid run("hiad", 30);
	//
	// utility tests
	run(test_sparse_cppad2eigen,            "sparse_cppad2eigen");
	run(test_solve_lower_cppad,             "solve_lower_cppad");
	run(test_jac_constraint,                "jac_constraint");
	run(test_control_jac_constraint,        "control::jac_contraint");
	//
	// control tests
	run(test_control_objective,             "control::objective");
	run(test_control_constraint,            "control::contraint");
	run(test_control_rec_constraint,        "control::rec_contraint");
	run(test_control_full_newton,           "control::full_newton");
	run(test_control_implicit_solver,       "control::implicit_solver");
	run(test_control_ad_rec_constraint,     "control::ad_rec_constraint");
	run(test_control_reduced_objective,     "control::reduced_objective");
	//
	// implicit kedem and newton tests
	run(test_cricle_implicit_kedem,         "circle_implicit_kedem");
	run(test_cricle_implicit_newton,        "circle_implicit_newton");
    //
    // check for memory leak
    bool memory_ok = CppAD::thread_alloc::free_all();
    // print summary at end
    bool ok = run.summary(memory_ok);
	// return
	if( ok )
		return 0;
	return 1;
}
