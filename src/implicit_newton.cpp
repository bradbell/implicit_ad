# include "implicit_newton.hpp"
# include "control.hpp"
/*
------------------------------------------------------------------------------
$begin test_circle_implicit_newton$$

$section Example / Test of Implicit Newton Class$$

$srcfile%src/implicit_newton.cpp%
	0%// BEGIN_TEST_CIRCLE_IMPLICIT_NEWTON%// END_TEST_CIRCLE_IMPLICIT_NEWTON%
1%$$

$end
*/
// BEGIN_TEST_CIRCLE_IMPLICIT_NEWTON
//
// simple implicit function y(x) defined by
// 0 = L(x, y) = x^2  + y(x)^2  - 1
//
class newton_circle {
public:
	newton_circle(void)
	{ }
	VECTOR(double) function(const VECTOR(double)& x)
	{	assert( x.size() == 1 );
		VECTOR(double) y(1);
		y[0] = sqrt( 1.0 - x[0] * x[0] );
		return y;
	}
	template <class Scalar>
	VECTOR(Scalar) derivative(const VECTOR(Scalar)& x, const VECTOR(Scalar)& y)
	{	// b = L_y (x, y)
		VECTOR(Scalar) b(1);
		b[0] = 2.0 * y[0];
		return b;
	}
	VECTOR( CppAD::AD<double> ) linear(
		const VECTOR( CppAD::AD<double>)&   b ,
		const VECTOR( CppAD::AD<double>)&   v )
	{	VECTOR( CppAD::AD<double> ) u(1);
		assert( v.size() == 1 );
		u[0] = v[0] / b[0];
		return u;
	}
};
bool test_cricle_implicit_newton(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	typedef CppAD::AD<double>   adouble;
	typedef CppAD::AD<adouble> a2double;

	// record z = L(x, y)
	VECTOR(a2double) a2xy(2), a2z(1);
	a2xy[0] = a2double( 0.5 );
	a2xy[1] = a2double( 0.5 );
	CppAD::Independent( a2xy );
	a2z[0] = a2xy[0] * a2xy[0] + a2xy[1] * a2xy[1] - 1.0;
	CppAD::ADFun<adouble> aL_fun(a2xy, a2z);

	// record F(x, y) = y
	VECTOR(adouble) axy(2), az(1);
	axy[0] = adouble( 0.5 );
	axy[1] = adouble( 0.5 );
	CppAD::Independent( axy );
	az[0] = axy[1];
	CppAD::ADFun<double> F_fun(axy, az);

	// solver used by newton_ad object
	newton_circle solve;

	// create implicit function object
	bool  full_step = false;
	size_t num_step = 2;
	implicit_newton<newton_circle> newton_ad(
		full_step, num_step, aL_fun, F_fun, solve
	);

	// Taylor coefficient vectors
	VECTOR(double) xk(1), yk(1);

	// zero order, y(t) = sqrt(1 - x(t) * x(t) )
	size_t k  = 0;
	double x0 = 0.5;
	xk[0]     = x0;
	yk        = newton_ad.Forward(k, xk);
	double y0 = sqrt(1.0 - x0 * x0);
	ok       &= CppAD::NearEqual(yk[0], y0, eps99, eps99);

	// differentiate zero order forward
	// y0 = sqrt( 1.0 - x0 * x0 )
	size_t q        = 1;
	VECTOR(double) w(q), dw(q);
	w[0]            = 1.0;
	dw              = newton_ad.Reverse(q, w);
	double dy0_dx0  = - x0 / y0;
	ok             &= CppAD::NearEqual(dw[0], dy0_dx0, eps99, eps99);

	// first order, y'(t) = - x(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	k         = 1;
	double x1 = 0.75;
	xk[0]     = x1;
	yk        = newton_ad.Forward(k, xk);
	double y1 = - x0 * x1 / y0;
	ok       &= CppAD::NearEqual(yk[0], y1, eps99, eps99);

	// differentiate first order forward
	// y1(x0, x1) = - x0 * x1 / y0
	q               = 2;
	w.resize(q);
	dw.resize(q);
	w[0]            = 0.0;
	w[1]            = 1.0;
	dw              = newton_ad.Reverse(q, w);
	double dy1_dx1  = - x0 / y0;
	ok           &= CppAD::NearEqual(dw[1], dy1_dx1, eps99, eps99);
	double dy1_dx0  = - x1 / y0 + ( x0 * x1 / (y0 * y0) ) * dy0_dx0;
	ok           &= CppAD::NearEqual(dw[0], dy1_dx0, eps99, eps99);

	// second order
	// y''(t) = - x'(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	//          - x(t) * x''(t) / sqrt(r0 - x(t) * x(t) )
	//          - (x(t) * x'(t))^2  / sqrt(r0 - x(t) * x(t) )^3
	k         = 2;
	double x2 = 0.25;
	xk[0]     = x2;
	yk        = newton_ad.Forward(k, xk);
	double y2 = - x1 * x1 / (2.0 * y0);
	y2       -= x0 * x2 / y0;
	y2       -= (x0 * x1) * (x0 * x1) / (2.0 * y0 * y0 * y0);
	ok       &= CppAD::NearEqual(yk[0], y2, eps99, eps99);
	//
	return ok;
}
// END_TEST_CIRCLE_IMPLICIT_NEWTON
