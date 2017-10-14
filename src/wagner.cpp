/*
$begin inverse$$
$spell
	Taylor
$$
$latex \newcommand{\B}[1]{{\bf #1}}$$

$section Wagner Inverse Function Derivative Algorithm$$

$head Notation$$
We can compute a smooth function
$latex Y : \B{R}^n \rightarrow \B{R}^n$$
and define the inverse function
$latex Z : \B{R}^n \rightarrow \B{R}^n$$; i.e.,
$latex Y[ Z(x) ] = x$$.
In an abuse of notation,
we fix a function $latex x : \B{R} \rightarrow \B{R}^n$$
and define the function $latex y : \B{R} \rightarrow \B{R}^n$$ by
$latex \[
	y(t) = Y( Z[ x(t) ] ] = x(t)
\] $$

$head Problem$$
We are given the Taylor coefficients
$latex \[
	x^k = x^{(k)} (0) / k!
\] $$
for $latex k = 0 , \ldots $$.
Our problem is to compute the coefficients
$latex \[
	z^k = z^{(k)} (0) / k!
\] $$
for $latex k = 0 , \ldots $$.

$head Zero Order$$
The value $latex z^0$$ is found by solving the equation
$latex Y( z^0 ) = x^0$$.

$head First Order$$
The value $latex z^1$$ is found by solving the equation
$latex \[
	Y^{(1)} ( z^0 ) z^1 = x^1
\] $$
In other words
$latex \[
	z^1 = Y^{(1)} ( z^0 )^{-1} x^1
\] $$

$head Induction$$
Suppose we know $latex z^j$$ for $latex z = 0 , \ldots , k-1$$
where $latex k \geq 1$$.
There is a function $latex H_k$$ such that
for any value of $latex z^k$$,
$latex \[
x^k
=
H_k ( z^0 , \ldots , z^{k-1} )
+
Y^{(1)} ( z^0 ) z^k
\]$$
We can use forward mode to compute
the value of $latex x^k$$ that would correspond to $latex z^k = 0$$; i.e.,
$latex \[
	H_k ( z^0 , \ldots , z^{k-1} )
\] $$
We can then solve for
$latex \[
z^k = Y^{(1)} ( z^0 )^{-1}
	[ x^k - H_k ( z^0 , \ldots , z^{k-1} ) ]
\] $$

$head Example$$
$srcfile%wagner.cpp%0%// BEGIN INVERSE%// END INVERSE%1%$$

$end
*/
// BEGIN INVERSE
# include <cppad/cppad.hpp>

bool inverse(void)
{	bool ok = true;
	//
	using CppAD::AD;
	using CppAD::vector;
	using CppAD::NearEqual;
	//
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// Case where Z(x) = asin(x) and Y(z) = sin(z)
	size_t n = 1, m = 1;
	vector< AD<double> > az(n), ay(m);
	az[0] = 1.0;
	CppAD::Independent( az );
	ay[0] = sin(az[0]);
	CppAD::ADFun<double> Y(az, ay);

	// zero order coefficients for x and z
	double x0 = 0.5;
	double z0 = asin(0.5);

	// zero order forward at x0
	vector<double> z(n), y(m);
	z[0] = z0;
	y    = Y.Forward(0, z);
	// should be near x0
	ok &= NearEqual(y[0], x0, eps99, eps99);

	// evaluate Y_z = Y^{(1)} ( z_0 )
	z[0] = 1.0;
	y    = Y.Forward(1, z);
	double Y_z = y[0];

	// solve for z1 in terms of x1
	double x1  = 0.25;
	double z1  = (1.0 / Y_z) * x1;

	// check z'(t) = z1
	// z'(t) =  x'(t) / sqrt( 1 - x(t) * x(t) )
	double root = sqrt(1.0 - x0 * x0);
	double zp   = x1 / root;
	ok &= NearEqual(z1, zp, eps99, eps99);

	// first order forward mode
	z[0] = z1;
	y    = Y.Forward(1, z);
	// should be near z1
	ok   &= NearEqual(y[0], x1, eps99, eps99);

	// second order forward with z2 = 0
	z[0] = 0.0;
	y    = Y.Forward(2, z);
	double H_2 = y[0];

	// solve for proper value of z2
	double x2  = 0.30;
	double z2  = (1.0 / Y_z) * (x2 - H_2);

	// check z''(t) / 2! = z2
	// z''(t) =  x''(t) / sqrt( 1 - x(t) * x(t) )
	//        + x'(t) * x'(t) * x(t) ( 1 - x(t) * x(t) )^{-3/2}
	double zpp = 2.0 * x2 / root;
	zpp       += x1 * x1 * x0 / (root * root * root);
	ok        &= NearEqual(zpp / 2.0, z2, eps99, eps99);

	// second order forward mode
	z[0] = z2;
	y    = Y.Forward(2, z);
	// should be near x2
	ok   &= NearEqual(y[0], x2, eps99, eps99);

	return ok;
}
// END INVERSE


/*
$begin implicit$$
$spell
	Taylor
$$
$latex \newcommand{\B}[1]{{\bf #1}}$$

$section Implicit Function Derivatives: Algorithm Similar to Wagner 2010$$

$head Notation$$
We are given a smooth function
$latex L : \B{R}^n \times \B{R}^m \rightarrow \B{R}^m $$
and define the implicit function
$latex Y : \B{R}^n \rightarrow \B{R}^m $$
by $latex L[ x , Y(x) ] = 0 $$.
In an abuse of notation,
we fix a function $latex x : \B{R} \rightarrow \B{R}^n$$
and define the function $latex y : \B{R} \rightarrow \B{R}^m$$ by
$latex \[
	L[ x(t) , y(t) ] = 0
\] $$

$head Problem$$
We are given the Taylor coefficients
$latex \[
	x^k = x^{(k)} (0) / k!
\] $$
for $latex k = 0 , \ldots $$.
Our problem is to compute the coefficients
$latex \[
	y^k = y^{(k)} (0) / k!
\] $$
for $latex k = 0 , \ldots $$.

$head Zero Order$$
The value $latex y^0$$ is found by solving the equation
$latex L( x^0 , y^0 ) = 0$$.

$head First Order$$
The value $latex y^1$$ is found by solving the equation
$latex \[
	L_x ( x^0 , y^0 ) x^1 + L_y ( x^0 , y^0 ) y^1 = 0
\] $$
In other words
$latex \[
	y^1 = - L_y ( x^0 , y^0 )^{-1} L_x ( x^0 , y^0 ) x^1
\] $$

$head Induction$$
Suppose we know $latex y^j$$ for $latex j = 0 , \ldots , k-1$$
where $latex k \geq 1$$.
There is a function $latex H_k$$ such that
for any value of $latex y^k$$,
$latex \[
v^k
=
H_k ( x^0 , \ldots , x^k , y^0 , \ldots , y^{k-1} )
+
L_y ( x^0 , y^0 ) y^k
\]$$
We can use forward mode to compute
the $latex v^k$$ that would correspond to $latex y^k = 0$$; i.e.,
$latex \[
	H_k ( x^0 , \ldots , x^k , y^0 , \ldots , y^{k-1} )
\] $$
We can then solve for
$latex \[
y^k = - L_y ( x^0 , y^0 )^{-1}
	H_k ( x^0 , \ldots , x^k , y^0 , \ldots , y^{k-1} )
\] $$

$head Example$$
$srcfile%wagner.cpp%0%// BEGIN IMPLICIT%// END IMPLICIT%1%$$

$end
*/
// BEGIN IMPLICIT
# include <cppad/cppad.hpp>

bool implicit(void)
{	bool ok = true;
	//
	using CppAD::AD;
	using CppAD::vector;
	using CppAD::NearEqual;
	//
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// L(x, y) = x^2 + y^2 = 1.0
	size_t n = 2, m = 1;
	vector< AD<double> > au(n), av(m);
	au[0] = 1.0;
	au[1] = 1.0;
	CppAD::Independent( au );
	AD<double> ax = au[0];
	AD<double> ay = au[1];
	av[0]         = 1.0 - (ax * ax + ay * ay);
	CppAD::ADFun<double> L(au, av);

	// zero order coefficients for x and y
	// y(t) = sqrt( 1 - x(t) * x(t) )
	double x0 = 0.5;
	double y0 = std::sqrt( 1.0 - x0 * x0 );

	// zero order forward at (x0, y0)
	vector<double> u(n), v(m);
	u[0] = x0;
	u[1] = y0;
	v    = L.Forward(0, u);
	// should be near zero
	ok &= NearEqual(v[0], 0.0, eps99, eps99);

	// evaluate L_y
	u[0] = 0.0;
	u[1] = 1.0;
	v    = L.Forward(1, u);
	double L_y = v[0];

	// first order forward mode with
	// x(t) = x0 + x1 * t, y(t) = y0
	double x1  = 0.5;
	u[0]       = x1;
	u[1]       = 0.0;
	v          = L.Forward(1, u);
	double H_1 = v[0];

	// solve for y1
	double y1 = - (1.0 / L_y) * H_1;

	// check y'(t) = y1
	// y'(t) = - x(t) * x'(t) / sqrt( 1 - x(t) * x(t) )
	double yp = - x0 * x1 / y0;
	ok &= NearEqual(y1, yp, eps99, eps99);

	// redo first order forward mode with proper value for y1
	u[0]  = x1;
	u[1]  = y1;
	v     = L.Forward(1, u);
	// should be near zero
	ok   &= NearEqual(v[0], 0.0, eps99, eps99);

	// second order forward mode with
	// x(t) = x0 + x1 * t + x2 * t^2, y(t) = y0 + y1 * t
	double x2  = 0.30;
	u[0]       = x2;
	u[1]       = 0.0;
	v          = L.Forward(2, u);
	double H_2 = v[0];

	// solve for y2
	double y2 = - (1.0 / L_y) * H_2;

	// check y''(t) / 2! = y2
	// y''(t) = - x'(t) * x'(t) / sqrt( 1 - x(t) * x(t) )
	//         - x(t) * x''(t) / sqrt( 1 - x(t) * x(t) )
	//         - (x(t) * x'(t))^2 / sqrt(1 - x(t) * x(t))^3
	double ypp = - x0 * x1 / y0 - x0 * 2.0 * x2 / y0;
	ypp       -= (x0 * x1) * (x0 * x1) / (y0 * y0 * y0);
	//
	ok   &= NearEqual(y2, ypp / 2.0, eps99, eps99);

	// redo second order forward mode with proper value for y2
	u[0]  = x2;
	u[1]  = y2;
	v     = L.Forward(2, u);
	// should be near zero
	ok   &= NearEqual(v[0], 0.0, eps99, eps99);
	//
	return ok;
}
// END IMPLICIT

int main(void)
{	bool ok = true;
	ok     &= inverse();
	ok     &= implicit();
	//
	if( ok )
	{	std::cout << "wagner: OK\n";
		return 0;
	}
	//
	std::cout << "wagner: Error\n";
	return 1;
}
