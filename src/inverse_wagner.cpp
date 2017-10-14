# include "inverse_wagner.hpp"

namespace { // BEGIN_EMPTY_NAMESPACE
//
typedef typename Eigen::Matrix<double, Eigen::Dynamic, 1>             vector;
typedef typename Eigen::Matrix< CppAD::AD<double>, Eigen::Dynamic, 1> a_vector;

// ----------------------------------------------------------------------------
// test a simple inversion of z = sin(x)
// ----------------------------------------------------------------------------
class solve_sin_x {
private:
	// f^{(1)} ( x^0 )
	double fp_;
public:
	solve_sin_x(void)
	{	fp_ = std::numeric_limits<double>::quiet_NaN();
	}
	vector function(const vector& z)
	{	vector x(1);
		assert( z.size() == 1 );
		x[0] = std::asin( z[0] );
		fp_  = std::cos( x[0] );
		return x;
	}
	vector derivative(const vector& v)
	{	vector u(1);
		assert( v.size() == 1 );
		u[0] = v[0] / fp_;
		return u;
	}
};
bool test_sin_x(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// record z = f(x)
	a_vector ax(1), az(1);
	ax[0] = 1.0;
	CppAD::Independent( ax );
	az[0] = sin( ax[0] );
	CppAD::ADFun<double> fun(ax, az);

	// solver object used by inverse function object
	solve_sin_x solve;

	// create inverse function object
	inverse_wagner<solve_sin_x> inv_fun(fun, solve);

	// Taylor coefficient vectors
	vector xk(1), zk(1);

	// zero order, x(t) = asin( z(t) )
	// x0 = x(0)
	size_t k  = 0;
	double z0 = 0.5;
	zk[0]     = z0;
	xk        = inv_fun.forward(k, zk);
	double x0 = asin(z0);
	ok       &= CppAD::NearEqual(xk[0], x0, eps99, eps99);

	// used to simplify expressions below
	double r0 = sqrt( 1 - z0 * z0 );


	// first order, x'(t) = z'(t) / sqrt( 1 - z(t) * z(t) )
	// x1 = x'(0)
	k         = 1;
	double z1 = 0.25;
	zk[0]     = z1;
	xk        = inv_fun.forward(k, zk);
	double x1 = z1 / r0;
	ok       &= CppAD::NearEqual(xk[0], x1, eps99, eps99);

	// second order,
	// x''(t) = z''(t) / sqrt( 1 - z(t) * z(t) )
	//        + z(t) * z'(t) * z'(t) / sqrt( 1 - z(t) * z(t )^3
	// x2     = x''(0) / 2
	k         = 2;
	double z2 = 0.75;
	zk[0]     = z2;
	xk        = inv_fun.forward(k, zk);
	double x2 = z2 / r0 + z0 * z1 * z1 / (2.0 * r0 * r0 * r0);
	ok       &= CppAD::NearEqual(xk[0], x2, eps99, eps99);

	//
	return ok;
}
// ----------------------------------------------------------------------------
// simple implicit function y(x) defined by x^2  + y(x)^2  = 1.0
// z = f(x, y) = (x, x^2 - y^2)
// f^{-1} (x, r) = [ x, sqrt(r - x * x) ]
// f^{-1} (x, 1) = [ x, y(x) ]
// ----------------------------------------------------------------------------
class solve_circle {
private:
	typedef typename
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
	// f^{(1)} ( x^0, y^0 )
	matrix inv_fp_;
public:
	solve_circle(void) : inv_fp_(2, 2)
	{	for(size_t i = 0; i < 2; i++)
			for(size_t j = 0; j < 2; j++)
				inv_fp_(i, j) = std::numeric_limits<double>::quiet_NaN();
	}
	vector function(const vector& z)
	{	vector u(2);
		assert( z.size() == 2 );
		double x0 = z[0];
		double r0 = z[1];
		double y0 = std::sqrt(r0 - x0 * x0);
		matrix fp(2, 2);
		// partial x w.r.t x
		fp(0, 0) = 1.0;
		// partial x w.r.t y
		fp(0, 1) = 0.0;
		// partial x^2 + y^2 w.r.t x
		fp(1, 0) = 2.0 * x0;
		// partial x^2 + y^2 w.r.t y
		fp(1, 1) = 2.0 * y0;
		//
		inv_fp_ = fp.inverse();
		//
		u[0] = x0;
		u[1] = y0;
		//
		return u;
	}
	vector derivative(const vector& v)
	{	vector u(2);
		assert( v.size() == 2 );
		u = inv_fp_ * v;
		return u;
	}
};
bool test_circle(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// record z = f(x, y)
	a_vector au(2), az(2);
	au[0] = 0.5;
	au[1] = 0.5;
	CppAD::Independent( au );
	az[0] = au[0];
	az[1] = au[0] * au[0] + au[1] * au[1];
	CppAD::ADFun<double> fun(au, az);

	// solver object used by inverse function object
	solve_circle solve;

	// create inverse function object
	inverse_wagner<solve_circle> inv_fun(fun, solve);

	// Taylor coefficient vectors
	vector uk(2), zk(2);

	// zero order, y(t) = sqrt(r0 - x(t) * x(t) )
	size_t k  = 0;
	double x0 = 0.5;
	double r0 = 1.0;
	zk[0]     = x0;
	zk[1]     = r0;
	uk        = inv_fun.forward(k, zk);
	double y0 = sqrt(r0 - x0 * x0);
	ok       &= CppAD::NearEqual(uk[1], y0, eps99, eps99);

	// first order, y'(t) = - x(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	k         = 1;
	double x1 = 0.75;
	zk[0]     = x1;
	zk[1]     = 0.0;
	uk        = inv_fun.forward(k, zk);
	double y1 = - x0 * x1 / y0;
	ok       &= CppAD::NearEqual(uk[1], y1, eps99, eps99);

	// second order
	// y''(t) = - x'(t) * x'(t) / sqrt(r0 - x(t) * x(t) )
	//          - x(t) * x''(t) / sqrt(r0 - x(t) * x(t) )
	//          - (x(t) * x'(t))^2  / sqrt(r0 - x(t) * x(t) )^3
	k         = 2;
	double x2 = 0.25;
	zk[0]     = x2;
	zk[1]     = 0.0;
	uk        = inv_fun.forward(k, zk);
	double y2 = - x1 * x1 / (2.0 * y0);
	y2       -= x0 * x2 / y0;
	y2       -= (x0 * x1) * (x0 * x1) / (2.0 * y0 * y0 * y0);
	ok       &= CppAD::NearEqual(uk[1], y2, eps99, eps99);
	//
	return ok;
}
} // END_EMPTY_NAMESPACE
// ---------------------------------------------------------------------------
// run tests
// ---------------------------------------------------------------------------
int main(void)
{	bool ok = true;
	ok     &= test_sin_x();
	ok     &= test_circle();
	//
	if( ok )
	{	std::cout << "inverse_wagner: OK\n";
		return 0;
	}
	//
	std::cout << "inverse_wagner: Error\n";
	return 1;
}
