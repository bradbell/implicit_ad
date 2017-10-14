# include "implicit_kedem.hpp"
# include "control.hpp"
/*
------------------------------------------------------------------------------
$begin test_circle_implicit_kedem$$

$section Example / Test of Implicit Wagner Class$$

$srcfile%implicit_kedem.cpp%
	0%// BEGIN_TEST_CIRCLE_IMPLICIT_WAGNER%// END_TEST_CIRCLE_IMPLICIT_WAGNER%
1%$$

$end
*/
// BEGIN_TEST_CIRCLE_IMPLICIT_WAGNER
//
// simple implicit function y(x) defined by
// 0 = L(x, y) = x^2  + y(x)^2  - 1
//
class kedem_circle {
public:
	kedem_circle(void)
	{ }
	VECTOR(double) function(const VECTOR(double)& x)
	{	assert( x.size() == 1 );
		VECTOR(double) y(1);
		y[0] = sqrt( 1.0 - x[0] * x[0] );
		return y;
	}
	VECTOR(double) derivative(const VECTOR(double)& x, const VECTOR(double)& y)
	{	// b = L_y (x, y)
		VECTOR(double) b(1);
		b[0] = 2.0 * y[0];
		return b;
	}
	VECTOR(double) linear(const VECTOR(double)& b, const VECTOR(double)& v)
	{	VECTOR(double) u(1);
		assert( v.size() == 1 );
		u[0] = v[0] / b[0];
		return u;
	}
};
bool test_cricle_implicit_kedem(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// record z = L(x, y)
	VECTOR( CppAD::AD<double> ) axy(2), az(1);
	axy[0] = 0.5;
	axy[1] = 0.5;
	CppAD::Independent( axy );
	az[0] = axy[0] * axy[0] + axy[1] * axy[1] - 1.0;
	CppAD::ADFun<double> L_fun(axy, az);

	// record F(x, y) = y
	CppAD::Independent( axy );
	az[0] = axy[1];
	CppAD::ADFun<double> F_fun(axy, az);

	// solver used by kedem_ad object
	kedem_circle solve;

	// create implicit function object
	implicit_kedem<kedem_circle> kedem_ad(L_fun, F_fun, solve);

	// Taylor coefficient vectors
	VECTOR(double) xk(1), yk(1);

	// zero order, y(t) = sqrt(1 - x(t) * x(t) )
	size_t k  = 0;
	double x0 = 0.5;
	xk[0]     = x0;
	bool store = true;
	yk         = kedem_ad.Forward(k, xk, store);
	double y0 = sqrt(1.0 - x0 * x0);
	ok       &= CppAD::NearEqual(yk[0], y0, eps99, eps99);

	// first order, y'(t) = - x(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	k         = 1;
	double x1 = 0.75;
	xk[0]     = x1;
	store     = true;
	yk        = kedem_ad.Forward(k, xk, store);
	double y1 = - x0 * x1 / y0;
	ok       &= CppAD::NearEqual(yk[0], y1, eps99, eps99);

	// second order
	// y''(t) = - x'(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	//          - x(t) * x''(t) / sqrt(r0 - x(t) * x(t) )
	//          - (x(t) * x'(t))^2  / sqrt(r0 - x(t) * x(t) )^3
	k         = 2;
	double x2 = 0.25;
	xk[0]     = x2;
	store     = false;
	yk        = kedem_ad.Forward(k, xk, false);
	double y2 = - x1 * x1 / (2.0 * y0);
	y2       -= x0 * x2 / y0;
	y2       -= (x0 * x1) * (x0 * x1) / (2.0 * y0 * y0 * y0);
	ok       &= CppAD::NearEqual(yk[0], y2, eps99, eps99);
	//
	return ok;
}
// END_TEST_CIRCLE_IMPLICIT_WAGNER
