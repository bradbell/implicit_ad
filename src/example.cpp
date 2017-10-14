# include <cppad/cppad.hpp>

namespace {
	using CppAD::AD;
	using CppAD::vector;

	// L(x, y)
	template <class Float>
	Float L(const Float& x, const Float& y)
	{	return Float(1.0) - x * x - y * y; }

	// L_y(x, y)
	template <class Float>
	Float L_y(const Float& x, const Float& y)
	{	return - Float(2.0) * y; }

	// half_step(z, x)
	template <class Float>
	Float half_step(double L_y_inv, const Float& x, const Float& z)
	{	return z - L_y_inv * L(x, z); }

	// full_step(z, x)
	template <class Float>
	Float full_step(const Float& x, const Float& z)
	{	return z - L(x, z) / L_y(x, z); }

	// Y(x) = ( 1 - x * x )^{-1/2}
	double Y0(double x)
	{	return CppAD::sqrt( 1.0 - x *  x ); }

	// Y'(x) = - x * (1 - x * x)^{-1/2}
	//       = - x / Y(x)
	double Y1(double x)
	{	double y  = Y0(x);
		return  - x / y;
	}

	// Y^{(2)} (x) = - [ Y(x) - x * Y'(x) ] / [ Y(x) * Y(x) ]
	//             = - Y(x)^{-1} - x * x * Y(x)^{-3}
	double Y2(double x)
	{	double y   = Y0(x);
		double ret = - ( 1.0 / y + x * x / (y * y * y) );
		ret       /=  2.0; // Y^{(2)} (x) / 2 !
		return ret;
	}

	// Y^{(3)} (x)
	// = Y(x)^{-2} * Y'(x) - 2 * x * Y(x)^{-3} + 3 * x * x * Y(x)^{-4} * Y'(x)
	// = - x * Y(x)^{-3} - 2 * x * Y(x)^{-3} - 3 * x * x * x * Y(x)^{-5}
	// = - 3 * x * Y(x)^{-3} [ 1 + x * x * Y(x)^{-2} ]
	double Y3(double x)
	{	double y   = Y0(x);
		double ret = - ( 3.0 * x * ( 1.0 + x * x / (y * y) ) / (y * y * y) );
		ret       /=  6.0; // Y^{(3)} (x) / 3 !
		return ret;
	}
}

int main(void)
{	bool ok    = true;
	double eps = 10. * std::numeric_limits<double>::epsilon();

	// vector and variables
	double x       = 0.25;
	double y       = Y0(x);
	double L_y_inv = 1.0 / L_y(x, y);

	// Record two full Newton steps as a function of x
	vector< AD<double> > ax(1), az(1);
	ax[0] = x;
	Independent(ax);
	az[0] = y;
	az[0] = full_step(ax[0], az[0]);
	az[0] = full_step(ax[0], az[0]);
	CppAD::ADFun<double> full(ax, az);

	// Record two half Newton steps as a function of x
	Independent(ax);
	az[0] = y;
	az[0] = half_step(L_y_inv, ax[0], az[0]);
	az[0] = half_step(L_y_inv, ax[0], az[0]);
	CppAD::ADFun<double> half(ax, az);

	// check zero order Taylor coefficient
	vector<double> x0(1), z0(1);
	double check;
	x0[0] = x;
	z0    = full.Forward(0, x0);
	ok &= fabs(1.0 - z0[0] / Y0(x) ) <= eps;
	z0    = half.Forward(0, x0);
	ok &= fabs(1.0 - z0[0] / Y0(x) ) <= eps;

	// check first order Taylor coefficient
	vector<double> x1(1), z1(1);
	x1[0] = 1.0;
	z1    = full.Forward(1, x1);
	ok &= fabs(1.0 - z1[0] / Y1(x) ) <= eps;
	z1    = half.Forward(1, x1);
	ok &= fabs(1.0 - z1[0] / Y1(x) ) <= eps;

	// check second order Taylor coefficient
	vector<double> x2(1), z2(1);
	x2[0]  = 0.0;
	z2     = full.Forward(2, x2);
	ok &= fabs(1.0 - z2[0] / Y2(x) ) <= eps;
	z2     = half.Forward(2, x2);
	ok &= fabs(1.0 - z2[0] / Y2(x) ) <= eps;

	// check third order Taylor coefficient
	vector<double> x3(1), z3(1);
	x3[0]  = 0.0;
	z3     = full.Forward(3, x3);
	ok &= fabs(1.0 - z3[0] / Y3(x) ) <= eps;

	// Note that two half steps does not give proper third derivatives
	z3     = half.Forward(3, x3);
	ok &= fabs(1.0 - z3[0] / Y3(x) ) > eps;

	if( ok )
	{	std::cout << "example: OK" << std::endl;
		return 0;
	}
	std::cout << "example: Error" << std::endl;
	return 1;
}


