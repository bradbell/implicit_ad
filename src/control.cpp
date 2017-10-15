# include "control.hpp"
# include "implicit_kedem.hpp"
# include "implicit_newton.hpp"
// ----------------------------------------------------------------------------
// BEGIN_TEST_OBJECTIVE
bool test_control_objective(void)
{	bool ok = true;
	double eps99   = 99.0 * std::numeric_limits<double>::epsilon();
	//
	VECTOR(double) q(4);
	MATRIX(double) x(2, 3);
	MATRIX(double) y(4, 3);
	for(size_t i = 0; i < 4; i++)
	{	q[i] = double(i + 1);
		for(size_t j = 0; j < 3; j++)
		{	y(i, j) = double(i + j + 1);
			if( i < 2 )
				x(i, j) = double(j - i + 2);
		}
	}
	double delta_t = 0.5;
	double F = control::objective(delta_t, q, x, y);
	//
	double check =  x(0,0) * x(0,0) + x(1,0) * x(1,0);
	check       +=  2.0 * ( x(0,1) * x(0,1) + x(1,1) * x(1,1) );
	check       +=  x(0,2) * x(0,2) + x(1,2) * x(1,2);
	check       *= delta_t / 4.0;
	for(size_t i = 0; i < 4; i++)
		check += q[i] * y(i,2) * y(i, 2) / 2.0;
	ok &= CppAD::NearEqual(F, check, eps99, eps99);
	//
	return ok;
}
// END_TEST_OBJECTIVE

// ----------------------------------------------------------------------------
// BEGIN_TEST_CONSTRAINT
bool test_control_constraint(void)
{	bool ok = true;
	double eps99   = 99.0 * std::numeric_limits<double>::epsilon();
	//
	VECTOR(double) p(4);
	MATRIX(double) x(2, 3);
	MATRIX(double) y(4, 3);
	for(size_t i = 0; i < 4; i++)
	{	p[i] = double(5 - i);
		for(size_t j = 0; j < 3; j++)
		{	y(i, j) = double(i + j + 1);
			if( i < 2 )
				x(i, j) = double(j - i + 2);
		}
	}
	double delta_t = 0.5;
	MATRIX(double) L = control::constraint(delta_t, p, x, y);
	ok &= L.rows() == 4;
	ok &= L.cols() == 3;
	//
	for(size_t i = 0; i < 4; i++)
	{	double check = y(i, 0) - p[i];
		ok &= CppAD::NearEqual(L(i, 0), check, eps99, eps99);
	}
	for(size_t j = 1; j < 3; j++)
	{	// use notation in Park2006 reference
		double u1 = x(0,j-1);
		double u2 = x(1,j-1);
		double x1 = y(0,j-1);
		double x2 = y(1,j-1);
		double x3 = y(2,j-1);
		double x4 = y(3,j-1);
		double r  = sqrt( (x1 + 1.0) * (x1 + 1.0) + x2 * x2 );
		double a  = 1.0 / (r * r * r) - 1.0;
		//
		double dxdt[4];
		dxdt[0] = x3;
		dxdt[1] = x4;
		dxdt[2] = + 2.0 * x4 - (1.0 + x1) * a + u1;
		dxdt[3] = - 2.0 * x3 -         x2 * a + u2;
		//
		for(size_t i = 0; i < 4; i++)
		{	double check = y(i,j) - y(i,j-1) - dxdt[i] * delta_t;
			ok &= CppAD::NearEqual(L(i, j), check, eps99, eps99);
		}
	}
	//
	return ok;
}
// END_TEST_CONSTRAINT
// ----------------------------------------------------------------------------
// BEGIN_TEST_REC_CONSTRAINT
bool test_control_rec_constraint(void)
{	bool ok = true;
	double eps99   = 99.0 * std::numeric_limits<double>::epsilon();
	//
	CppAD::ADFun<double>       L_fun;
	size_t                     J = 4;
	double                     delta_t = 0.5;
	VECTOR(double) p(4);
	for(size_t i = 0; i < 4; i++)
		p[i] = double(5 - i);
	control::rec_constraint(L_fun, J, delta_t, p);
	//
	// argument values
	MATRIX(double) x(2, J), y(4, J);
	for(size_t i = 0; i < 4; i++)
	{	for(size_t j = 0; j < J; j++)
		{	y(i, j) = double(2 * J + 4 * j + i);
			if( i < 2 )
				x(i, j) = double(2 * j + i);
		}
	}
	// pack (x, y) into xy_vec
	VECTOR(double) xy_vec(6 * J);
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 2; i++)
			xy_vec[ j * 2 + i] = x(i, j);
		for(size_t i = 0; i < 4; i++)
			xy_vec[ 2 * J + j * 4 + i] = y(i, j);
	}
	// compute constraint function using L_fun
	VECTOR(double) L_vec = L_fun.Forward(0, xy_vec);
	ok &= size_t( L_vec.size() ) == 4 * J;
	//
	// unpack L_vec
	MATRIX(double) L(4, J);
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 4; i++)
			L(i, j) = L_vec[ j * 4 + i];
	}
	//
	for(size_t i = 0; i < 4; i++)
	{	double check = y(i, 0) - p[i];
		ok &= CppAD::NearEqual(L(i, 0), check, eps99, eps99);
	}
	for(size_t j = 1; j < J; j++)
	{	// use notation in Park2006 reference
		double u1 = x(0,j-1);
		double u2 = x(1,j-1);
		double x1 = y(0,j-1);
		double x2 = y(1,j-1);
		double x3 = y(2,j-1);
		double x4 = y(3,j-1);
		double r  = sqrt( (x1 + 1.0) * (x1 + 1.0) + x2 * x2 );
		double a  = 1.0 / (r * r * r) - 1.0;
		//
		double dxdt[4];
		dxdt[0] = x3;
		dxdt[1] = x4;
		dxdt[2] = + 2.0 * x4 - (1.0 + x1) * a + u1;
		dxdt[3] = - 2.0 * x3 -         x2 * a + u2;
		//
		for(size_t i = 0; i < 4; i++)
		{	double check = y(i,j) - y(i,j-1) - dxdt[i] * delta_t;
			ok &= CppAD::NearEqual(L(i, j), check, eps99, eps99);
		}
	}
	//
	return ok;
}
// END_TEST_REC_CONSTRAINT
// ---------------------------------------------------------------------------
// BEGIN_TEST_FULL_NEWTON
bool test_control_full_newton(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	//
	// record L_fun
	CppAD::ADFun<double>       L_fun;
	size_t                     J = 10;
	double                     delta_t = 0.1;
	VECTOR(double) p(4);
	for(size_t i = 0; i < 4; i++)
		p[i] = 0.2;
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m      = L_fun.Range();
	size_t n      = L_fun.Domain() - m;
	//
	// value of x and initial y
	VECTOR(double) xy_in(n + m), xy_out;
	for(size_t i = 0; i < n+m; i++)
		xy_in[i] = 0.0;
	//
	double criteria  = double(n + m) * eps99;
	size_t max_itr   = 10;
	CPPAD_SPARSE(double)    L_y;
	CppAD::sparse_jac_work  work;
	jac_constraint(L_y, L_fun, xy_in, work);
	//
	size_t num_itr = control::full_newton(
		xy_out, xy_in, L_fun, criteria, max_itr, L_y, work
	);
	// max sure it did not require maximum nunber of iterations
	ok &= num_itr < max_itr;
	// make sure xy has not gone to a huge value
	ok &= norm_squared( xy_out ) / (6 * J) < 1.0;
	// check the convergence criteria
	VECTOR(double) L = L_fun.Forward(0, xy_out);
	ok &= norm_squared( L ) < criteria;
	//
	return ok;
}
// END_TEST_FULL_NEWTON
// ---------------------------------------------------------------------------
// BEGIN_TEST_CONTROL_IMPLICIT_SOLVER
bool test_control_implicit_solver(void)
{	bool ok      = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	//
	// record L_fun
	CppAD::ADFun<double>       L_fun;
	size_t                     J = 5;
	double                     delta_t = 0.1;
	VECTOR(double) p(4);
	for(size_t i = 0; i < 4; i++)
		p[i] = 0.2;
	control::rec_constraint(L_fun, J, delta_t, p);
	//
	size_t m      = L_fun.Range();
	size_t n      = L_fun.Domain() - m;
	//
	// create implicit solver
	double criteria  = double(n + m) * eps99;
	CppAD::ADFun< CppAD::AD<double> > not_used;
	control::implicit_solver solver(L_fun, not_used, criteria);
	//
	// set x
	VECTOR(double) x(n);
	for(size_t i = 0; i < n; i++)
		x[i] = 0.0;
	//
	// solve L(x, y) == 0
	VECTOR(double) y = solver.function(x);
	ok &= size_t( y.size() ) == m;
	//
	// check L(x, y) is near zero
	VECTOR(double) xy(n + m);
	join_vector(xy, x, y);
	VECTOR(double) L  = L_fun.Forward(0, xy);
	ok &= norm_squared(L) < criteria;
	//
	// solve L_y (x, y) * u = v
	VECTOR(double) v(m);
	for(size_t i = 0; i < m; i++)
		v[i] = double(i - 5.0);
	VECTOR(double) b = solver.derivative(x, y);
	VECTOR(double) u = solver.linear(b, v);
	//
	// compute L_y (x, y)
	CPPAD_SPARSE(double) L_y;
	CppAD::sparse_jac_work work;
	jac_constraint(L_y, L_fun, xy, work);
	SparseMatrix<double> L_y_eigen;
	sparse_cppad2eigen(L_y, L_y_eigen);
	//
	// check L_y (x, y) * u = v
	VECTOR(double) d = L_y_eigen * u - v;
	ok &= norm_squared(d) < criteria;
	//
	return ok;
}
// END_TEST_CONTROL_IMPLICIT_SOLVER
// ---------------------------------------------------------------------------
// BEGIN_TEST_AD_REC_CONSTRAINT
bool test_control_ad_rec_constraint(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	//
	// ---------------------------------------------------------------------
	// create afun
	typedef CppAD::AD<double> adouble;
	CppAD::ADFun<adouble>      afun;
	size_t                     J = 3;
	adouble                    adelta_t = 0.5;
	VECTOR(adouble) ap(4);
	for(size_t i = 0; i < 4; i++)
		ap[i] = adouble(5 - i);
	control::rec_constraint(afun, J, adelta_t, ap);
	size_t m = afun.Range();
	size_t n = afun.Domain() - m;
	// ----------------------------------------------------------------------
	// Create L_fun using afun
	VECTOR(adouble) axy(n + m), aL(m);
	for(size_t i = 0; i < n + m; i++)
		axy[i] = adouble(0.0);
	CppAD::Independent( axy );
	aL = afun.Forward(0, axy);
	CppAD::ADFun<double> L_fun(axy, aL);
	// ----------------------------------------------------------------------
	// check L_fun
	// ----------------------------------------------------------------------
	// argument values
	MATRIX(double) x(2, J);
	MATRIX(double) y(4, J);
	for(size_t i = 0; i < 4; i++)
	{	for(size_t j = 0; j < J; j++)
		{	y(i, j) = double(2 * J + 4 * j + i);
			if( i < 2 )
				x(i, j) = double(2 * j + i);
		}
	}
	// pack (x, y) into xy_vec
	VECTOR(double) xy_vec(6 * J);
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 2; i++)
			xy_vec[ j * 2 + i] = x(i, j);
		for(size_t i = 0; i < 4; i++)
			xy_vec[ 2 * J + j * 4 + i] = y(i, j);
	}
	// compute constraint function using L_fun
	VECTOR(double) L_vec = L_fun.Forward(0, xy_vec);
	ok &= size_t( L_vec.size() ) == 4 * J;
	//
	// unpack L_vec
	MATRIX(double) L(4, J);
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 4; i++)
			L(i, j) = L_vec[ j * 4 + i];
	}
	//
	for(size_t i = 0; i < 4; i++)
	{	double check = y(i, 0) - Value(ap[i]);
		ok &= CppAD::NearEqual(L(i, 0), check, eps99, eps99);
	}
	for(size_t j = 1; j < J; j++)
	{	// use notation in Park2006 reference
		double u1 = x(0,j-1);
		double u2 = x(1,j-1);
		double x1 = y(0,j-1);
		double x2 = y(1,j-1);
		double x3 = y(2,j-1);
		double x4 = y(3,j-1);
		double r  = sqrt( (x1 + 1.0) * (x1 + 1.0) + x2 * x2 );
		double a  = 1.0 / (r * r * r) - 1.0;
		//
		double dxdt[4];
		dxdt[0] = x3;
		dxdt[1] = x4;
		dxdt[2] = + 2.0 * x4 - (1.0 + x1) * a + u1;
		dxdt[3] = - 2.0 * x3 -         x2 * a + u2;
		//
		for(size_t i = 0; i < 4; i++)
		{	double check = y(i,j) - y(i,j-1) - dxdt[i] * Value(adelta_t);
			ok &= CppAD::NearEqual(L(i, j), check, eps99, eps99);
		}
	}
	//
	return ok;
}
// END_TEST_AD_REC_CONSTRAINT
/*
-- ---------------------------------------------------------------------------
$begin test_control_reduced_objective$$
$spell
	Kedem
$$

$section Example / Test of Control Problem Reduced Objective$$

$head Purpose$$
The reduced objective is defined by
$latex \[
	R(x) = F[ x , Y(x) ]
\]$$
This example / test compares derivatives of $latex R(x)$$
computed by the Kedem, Full Newton, and Partial Newton methods.

$head Source$$
$srcfile%src/control.cpp%
	0%// BEGIN_TEST_REDUCED_OBJECTIVE%// END_TEST_REDUCED_OBJECTIVE%
1%$$

$end
*/
// BEGIN_TEST_REDUCED_OBJECTIVE
bool test_control_reduced_objective(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	typedef CppAD::AD<double> adouble;
	// ----------------------------------------------------------------------
	// control problem parameters
	size_t   J = 3;
	double   delta_t = 0.1;
	VECTOR(double) p(4), q(4);
	for(size_t i = 0; i < 4; i++)
	{	p[i] = double(5 - i);
		q[i] = double(i + 1);
	}
	// ----------------------------------------------------------------------
	// L_fun, F_fun
	CppAD::ADFun<double> F_fun, L_fun;
	control::rec_objective(F_fun, J,  delta_t, q);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	ok &= F_fun.Range() == 1;
	ok &= F_fun.Domain() == n + m;
	// -----------------------------------------------------------------------
	// aL_fun
	CppAD::ADFun<adouble> aL_fun;
	adouble adelta_t = delta_t;
	VECTOR(adouble) ap(4);
	for(size_t i = 0; i < 4; i++)
		ap[i] = p[i];
	control::rec_constraint(aL_fun, J, adelta_t, ap);
	// -----------------------------------------------------------------------
	// control_solve
	double criteria = double(n + m) * eps99;
	control::implicit_solver control_solve(L_fun, aL_fun, criteria);
	// -----------------------------------------------------------------------
	// R_kedem, R_full, R_partial
	implicit_kedem<control::implicit_solver> R_kedem(
		L_fun, F_fun, control_solve
	);
	bool full_step  = true;
	size_t num_step = 2;
	implicit_newton<control::implicit_solver> R_full(
		full_step, num_step, aL_fun, F_fun, control_solve
	);
	full_step  = false;
	implicit_newton<control::implicit_solver> R_partial(
		full_step, num_step, aL_fun, F_fun, control_solve
	);
	// ----------------------------------------------------------------------
	// zero order forward
	VECTOR(double) x0(n);
	for(size_t i = 0; i < n; i++)
		x0[i] = 0.1;
	bool store = true;
	VECTOR(double) kedem_R0     = R_kedem.Forward(0, x0, store);
	VECTOR(double) full_R0      = R_full.Forward(0, x0);
	VECTOR(double) partial_R0   = R_partial.Forward(0, x0);
	ok &= CppAD::NearEqual(kedem_R0[0], partial_R0[0], eps99, eps99);
	ok &= CppAD::NearEqual(kedem_R0[0], full_R0[0],    eps99, eps99);
	// -----------------------------------------------------------------------
	// compute and check gradient w.r.t. x
	VECTOR(double) w(1);
	w[0] = 1.0;
	VECTOR(double) partial_dR = R_partial.Reverse(1, w);
	VECTOR(double) full_dR    = R_partial.Reverse(1, w);
	VECTOR(double) x1(n), kedem_dR(1);
	for(size_t i = 0; i < n; i++)
		x1[i] = 0.0;
	for(size_t j = 0; j < n; j++)
	{	x1[j] = 1.0;
		store = false;
		kedem_dR  = R_kedem.Forward(1, x1, store);
		ok &= CppAD::NearEqual(kedem_dR[0], partial_dR[j], eps99, eps99);
		ok &= CppAD::NearEqual(kedem_dR[0], full_dR[j],    eps99, eps99);
		x1[j] = 0.0;
	}
	// -----------------------------------------------------------------------
	// compute Hessian w.r.t. x corresponding to R_partial or R_full
	MATRIX(double) H(n, n);
	VECTOR(double) d_partial(n), d_full(n);
	w.resize(2);
	w[0] = 0.0;
	w[1] = 1.0;
	for(size_t i = 0; i < n; i++)
		x1[i] = 0.0;
	for(size_t i = 0; i < n; i++)
	{	x1[i] = 1.0;
		//
		// first order forward
		R_partial.Forward(1, x1);
		R_full.Forward(1, x1);
		//
		// second order reverse
		d_partial = R_partial.Reverse(2, w);
		d_full    = R_full.Reverse(2, w);
		for(size_t j = 0; j < n; j++)
		{	ok &= CppAD::NearEqual(
				d_partial[2 * j + 0], d_full[2 * j + 0], eps99, eps99
			);
			ok &= CppAD::NearEqual(
				d_partial[2 * j + 1], d_full[2 * j + 1], eps99, eps99
			);
			H(i, j) = d_partial[2 * j + 0];
		}
		x1[i] = 0.0;
	}
	// -----------------------------------------------------------------------
	// compare H with Hessian corresponding to kedem
	VECTOR(double) x2(n), ddR(1);
	for(size_t i = 0; i < n; i++)
		x2[i] = 0.0;
	for(size_t i = 0; i < n; i++)
	{	// check i-th diagonal element
		x1[i] = 1.0;
		//
		store = true;
		R_kedem.Forward(1, x1, true);
		store = false;
		ddR = R_kedem.Forward(2, x2, false);
		ok &= CppAD::NearEqual(ddR[0], H(i, i) / 2.0, eps99, eps99);
		//
		// check off diagonal elements
		for(size_t j = 0; j < n; j++) if( j != i )
		{	x1[j] = 1.0;
			//
			store = true;
			R_kedem.Forward(1, x1, store);
			store = false;
			ddR = R_kedem.Forward(2, x2, store);
			double Hij = ( 2.0 * ddR[0] - H(i, i) - H(j, j) ) / 2.0;
			ok &= CppAD::NearEqual(Hij, H(i, j), eps99, eps99);
			//
			x1[j] = 0.0;
		}
		//
		x1[i] = 0.0;
	}
	return ok;
}
// END_TEST_REDUCED_OBJECTIVE
