# include "control.hpp"
# include "implicit_kedem.hpp"
# include "implicit_newton.hpp"
# include <cppad/speed/uniform_01.hpp>
/*
-------------------------------------------------------------------------------
$begin set_T_p_and_q$$

$section Set T, p, and q$$

$head T$$
This is the total time in the optimal control example.

$head p$$
This is the initial value for the state vector at time zero.

$head q$$
This is the multiplier, in the objective,
for each component of the final state vector squared.

$srccode%cpp% */
void set_T_p_and_q(double& T, VECTOR(double)& p, VECTOR(double) &q )
{	T    = 1.0;
	//
	p[0] = 0.2;
	p[1] = 0.2;
	p[2] = 0.1;
	p[3] = 0.1;
	//
	q[0] = 25.0;
	q[1] = 25.0;
	q[2] = 25.0;
	q[3] = 25.0;
	//
	return;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin repeat_kedem_gradient$$
$spell
	Kedem
$$

$section Repeated Computation of Control Problem Gradient Using Kedem Method$$

$head Syntax$$
$codei%repeat_kedem_gradient(%repeat%, %J%, %x%, %grad%)
%$$
$codei%repeat_kedem_gradient(%repeat%, %size%)
%$$

$head repeat$$
is the number of times to choose a new value for the controls.

$head J$$
is the number of time points in the discrete version of the control problem.

$head x$$
If $icode x$$ has size zero, it is not used.
Otherwise, it must have size $codei%2*%J%$$.
In this case it specifies the value of the control vector for the
first gradient computation.

$head grad$$
This vector must have size $codei%2*%J%$$.
The input value of its elements does not matter.
Upon return, it contains the value last gradient computed
by $code repeat_kedem_gradient$$.

$head size$$
This is the number time points in the problem; i.e.,
the value of $icode J$$ for this test.

$srcfile%src/time.cpp%
	0%// BEGIN_REPEAT_KEDEM_GRADIENT%// END_REPEAT_KEDEM_GRADIENT%
1%$$

$end
*/
// BEGIN_REPEAT_KEDEM_GRADIENT
void repeat_kedem_gradient(
	size_t repeat           ,
	size_t J                ,
	const VECTOR(double)& x ,
	VECTOR(double)& grad    )
{	assert( x.size() == 0 || x.size() == grad.size() );
	assert( size_t( grad.size() ) == 2 * J );
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	// ----------------------------------------------------------------------
	// control problem parameters
	double T;
	VECTOR(double) p(4), q(4);
	set_T_p_and_q(T, p, q);
	double delta_t = T / double(J - 1);
	// ----------------------------------------------------------------------
	// L_fun, F_fun
	CppAD::ADFun<double> F_fun, L_fun;
	control::rec_objective(F_fun, J,  delta_t, q);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	assert( F_fun.Range() == 1 );
	assert( F_fun.Domain() == n + m );
	// -----------------------------------------------------------------------
	// control_solve
	double criteria = double(n + m) * eps99;
	CppAD::ADFun< CppAD::AD<double> > not_used;
	control::implicit_solver control_solve(L_fun, not_used, criteria);
	// -----------------------------------------------------------------------
	// R_kedem
	implicit_kedem<control::implicit_solver> R_kedem(
		L_fun, F_fun, control_solve
	);
	// x0, x1, dR
	VECTOR(double) x0(n), x1(n), dR(1);
	for(size_t i = 0; i < n; i++)
		x1[i] = 0.0;
	if( x.size() == 0 )
		CppAD::uniform_01(n, x0);
	else
		x0 = x;
	for(size_t count = 0; count < repeat; count++)
	{
		// zero order forward
		bool store = true;
		R_kedem.Forward(0, x0, store);
		//
		// use first order forward to calculate gradient
		for(size_t j = 0; j < n; j++)
		{	x1[j] = 1.0;
			store = false;
			dR    = R_kedem.Forward(1, x1, store);
			grad[j] = dR[0];
			x1[j] = 0.0;
		}
		// pick a random value for next x0
		CppAD::uniform_01(n, x0);
	}
}
void time_kedem_gradient(size_t size, size_t repeat)
{	size_t J = size;
	VECTOR(double) grad(2 * J), x(0);
	repeat_kedem_gradient(repeat, J, x, grad);
}
// END_REPEAT_KEDEM_GRADIENT
/*
-------------------------------------------------------------------------------
$begin repeat_newton_gradient$$
$spell
	Kedem
$$

$section Repeated Computation of Control Problem Gradient Using Newton Method$$

$head Syntax$$
$codei%repeat_newton_gradient(%repeat%, %J%, %reverse%, %x%, %grad%)
%$$
$codei%repeat_forward_gradient(%repeat%, %size%)
%$$
$codei%repeat_reverse_gradient(%repeat%, %size%)
%$$

$head repeat$$
is the number of times to choose a new value for the controls.

$head J$$
is the number of time points in the discrete version of the control problem.

$head reverse$$
If this is true, reverse mode is used to compute the gradient.
Otherwise, forward mode is used.

$head x$$
If $icode x$$ has size zero, it is not used.
Otherwise, it must have size $codei%2*%J%$$.
In this case it specifies the value of the control vector for the
first gradient computation.

$head grad$$
This vector must have size $codei%2*%J%$$.
The input value of its elements does not matter.
Upon return, it contains the value last gradient computed
by $code repeat_kedem_gradient$$.

$head size$$
This is the number time points in the problem; i.e.,
the value of $icode J$$ for this test.

$srcfile%src/time.cpp%
	0%// BEGIN_REPEAT_NEWTON_GRADIENT%// END_REPEAT_NEWTON_GRADIENT%
1%$$

$end
*/
// BEGIN_REPEAT_NEWTON_GRADIENT
void repeat_newton_gradient(
	size_t                repeat    ,
	size_t                J         ,
	bool                  reverse   ,
	const VECTOR(double)& x         ,
	VECTOR(double)&       grad      )
{	assert( x.size() == 0 || x.size() == grad.size() );
	assert( size_t( grad.size() ) == 2 * J );
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	typedef CppAD::AD<double> adouble;
	// ----------------------------------------------------------------------
	// control problem parameters
	double T;
	VECTOR(double) p(4), q(4);
	set_T_p_and_q(T, p, q);
	double delta_t = T / double(J - 1);
	// ----------------------------------------------------------------------
	// L_fun, F_fun
	CppAD::ADFun<double> F_fun, L_fun;
	control::rec_objective(F_fun, J,  delta_t, q);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	assert( F_fun.Range() == 1 );
	assert( F_fun.Domain() == n + m );
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
	// R_newton
	size_t num_step  = 1;
	bool   full_step = false; // does not matter when num_step = 1
	implicit_newton<control::implicit_solver> R_newton(
		full_step, num_step, aL_fun, F_fun, control_solve
	);
	// x0, w
	VECTOR(double) x0(n), w(1), x1(n), dR(1);
	w[0] = 1.0;
	for(size_t j = 0; j < n; j++)
		x1[j] = 0.0;
	if( x.size() == 0 )
		CppAD::uniform_01(n, x0);
	else
		x0 = x;
	for(size_t count = 0; count < repeat; count++)
	{
		//
		// zero order forward
		R_newton.Forward(0, x0);
		//
		if( reverse )
			grad = R_newton.Reverse(1, w);
		else
		{	// use first order forward to calculate gradient
			for(size_t j = 0; j < n; j++)
			{	x1[j] = 1.0;
				dR    = R_newton.Forward(1, x1);
				grad[j] = dR[0];
				x1[j] = 0.0;
			}
		}
		//
		// pick a random value for next x0
		CppAD::uniform_01(n, x0);
	}
}
void time_forward_gradient(size_t size, size_t repeat)
{	size_t J       = size;
	bool   reverse = false;
	VECTOR(double) grad(2 * J), x(0);
	repeat_newton_gradient(repeat, J, reverse, x, grad);
}
void time_reverse_gradient(size_t size, size_t repeat)
{	size_t J       = size;
	bool   reverse = true;
	VECTOR(double) grad(2 * J), x(0);
	repeat_newton_gradient(repeat, J, reverse, x, grad);
}
// END_REPEAT_NEWTON_GRADIENT

/*
-------------------------------------------------------------------------------
$begin repeat_kedem_hessian$$
$spell
	Kedem
	hess
$$

$section Repeated Computation of Control Problem Hessian Using Kedem Method$$

$head Syntax$$
$codei%repeat_kedem_hessian(%repeat%, %J%, %x%, %hess%)
%$$
$codei%repeat_kedem_hessian(%repeat%, %size%)
%$$

$head repeat$$
is the number of times to choose a new value for the controls.

$head J$$
is the number of time points in the discrete version of the control problem.

$head x$$
If $icode x$$ has size zero, it is not used.
Otherwise, it must have size $codei%2*%J%$$.
In this case it specifies the value of the control vector for the
first Hessian computation.

$head hess$$
This vector must have size $codei%4*%J%*%J%$$.
The input value of its elements does not matter.
Upon return, it contains the value last Hessian computed
by $code repeat_kedem_hessian$$.

$head size$$
This is the number time points in the problem; i.e.,
the value of $icode J$$ for this test.

$srcfile%src/time.cpp%
	0%// BEGIN_REPEAT_KEDEM_HESSIAN%// END_REPEAT_KEDEM_HESSIAN%
1%$$

$end
*/
// BEGIN_REPEAT_KEDEM_HESSIAN
void repeat_kedem_hessian(
	size_t                repeat  ,
	size_t                J       ,
	const VECTOR(double)& x       ,
	VECTOR(double)&       hess    )
{	assert( x.size() == 0 || size_t( x.size() ) == 2 * J );
	assert( size_t( hess.size() ) == 4 * J * J );
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	// ----------------------------------------------------------------------
	// control problem parameters
	double T;
	VECTOR(double) p(4), q(4);
	set_T_p_and_q(T, p, q);
	double delta_t = T / double(J - 1);
	// ----------------------------------------------------------------------
	// L_fun, F_fun
	CppAD::ADFun<double> F_fun, L_fun;
	control::rec_objective(F_fun, J,  delta_t, q);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	assert( F_fun.Range() == 1 );
	assert( F_fun.Domain() == n + m );
	// -----------------------------------------------------------------------
	// control_solve
	double criteria = double(n + m) * eps99;
	CppAD::ADFun< CppAD::AD<double> > not_used;
	control::implicit_solver control_solve(L_fun, not_used, criteria);
	// -----------------------------------------------------------------------
	// R_kedem
	implicit_kedem<control::implicit_solver> R_kedem(
		L_fun, F_fun, control_solve
	);
	//
	VECTOR(double) x0(n), x1(n), x2(n), ddR(1);
	for(size_t i = 0; i < n; i++)
	{	x1[i] = 0.0;
		x2[i] = 0.0;
	}
	if( x.size() == 0 )
		CppAD::uniform_01(n, x0);
	else
		x0 = x;
	for(size_t count = 0; count < repeat; count++)
	{
		// zero order forward
		bool store = true;
		R_kedem.Forward(0, x0, store);
		//
		// calculate the diagonal of Hessian
		for(size_t j = 0; j < n; j++)
		{	x1[j] = 1.0;
			store = true;
			R_kedem.Forward(1, x1, store);
			store = false;
			ddR             = R_kedem.Forward(2, x2, store);
			hess[j * n + j] = 2.0 * ddR[0];
			x1[j]           = 0.0;
		}
		//
		// calculate off diagonal elements
		for(size_t i = 0; i < n; i++)
		{	x1[i] = 1.0;
			for(size_t j = 0; j < i; j++) if( j != i )
			{	x1[j] = 1.0;
				//
				store = true;
				R_kedem.Forward(1, x1, store);
				store = false;
				ddR = R_kedem.Forward(2, x2, store);
				double Hii = hess[i * n + i];
				double Hjj = hess[j * n + j];
				double Hij = (2.0 * ddR[0] - Hii - Hjj) / 2.0;
				hess[i * n + j] = Hij;
				hess[j * n + i] = Hij;
				//
				x1[j]           = 0.0;
			}
			x1[i] = 0.0;
		}
		// pick a random value for nexgt x0
		CppAD::uniform_01(n, x0);
	}
}
void time_kedem_hessian(size_t size, size_t repeat)
{	size_t J = size;
	VECTOR(double) hess(2 * J * 2 * J), x(0);
	repeat_kedem_hessian(repeat, J, x, hess);
}
// END_REPEAT_KEDEM_HESSIAN
/*
-------------------------------------------------------------------------------
$begin repeat_newton_hessian$$
$spell
	Kedem
	hess
$$

$section Repeated Computation of Control Problem Hessian Using Newton Method$$

$head Syntax$$
$codei%repeat_kedem_hessian(%repeat%, %J%, %full_step%, %x%, %hess%)
%$$
$codei%repeat_partial_hessian(%repeat%, %size%)
%$$
$codei%repeat_full_hessian(%repeat%, %size%)
%$$

$head repeat$$
is the number of times to choose a new value for the controls.

$head J$$
is the number of time points in the discrete version of the control problem.

$head full_step$$
If this is true, the full Newton step is used to compute the Hessian.
Otherwise, the partial Newton step mode is used.

$head x$$
If $icode x$$ has size zero, it is not used.
Otherwise, it must have size $codei%2*%J%$$.
In this case it specifies the value of the control vector for the
first Hessian computation.

$head hess$$
This vector must have size $codei%4*%J%*%J%$$.
The input value of its elements does not matter.
Upon return, it contains the value last Hessian computed
by $code repeat_newton_hessian$$.

$head size$$
This is the number time points in the problem; i.e.,
the value of $icode J$$ for this test.

$srcfile%src/time.cpp%
	0%// BEGIN_REPEAT_NEWTON_HESSIAN%// END_REPEAT_NEWTON_HESSIAN%
1%$$

$end
*/
// BEGIN_REPEAT_NEWTON_HESSIAN
void repeat_newton_hessian(
	size_t                repeat     ,
	size_t                J          ,
	bool                  full_step	 ,
	const VECTOR(double)& x          ,
	VECTOR(double)&       hess       )
{	assert( x.size() == 0 || size_t( x.size() ) == 2 * J );
	assert( size_t( hess.size() ) == 2 * J * 2 * J );
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	typedef CppAD::AD<double> adouble;
	// ----------------------------------------------------------------------
	// control problem parameters
	double T;
	VECTOR(double) p(4), q(4);
	set_T_p_and_q(T, p, q);
	double delta_t = T / double(J - 1);
	// ----------------------------------------------------------------------
	// L_fun, F_fun
	CppAD::ADFun<double> F_fun, L_fun;
	control::rec_objective(F_fun, J,  delta_t, q);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	assert( F_fun.Range() == 1 );
	assert( F_fun.Domain() == n + m );
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
	// R_newton
	size_t num_step = 2;
	implicit_newton<control::implicit_solver> R_newton(
		full_step, num_step, aL_fun, F_fun, control_solve
	);
	//
	VECTOR(double) x0(n), x1(n), w(2), dRR(2 * n);
	w[0] = 0.0;
	w[1] = 1.0;
	for(size_t i = 0; i < n; i++)
		x1[i] = 0.0;
	if( x.size() == 0 )
		CppAD::uniform_01(n, x0);
	else
		x0 = x;
	for(size_t count = 0; count < repeat; count++)
	{	// zero order forward
		R_newton.Forward(0, x0);
		for(size_t i = 0; i < n; i++)
		{	x1[i] = 1.0;
			R_newton.Forward(1, x1);
			dRR = R_newton.Reverse(2, w);
			for(size_t j = 0; j < n; j++)
				hess[i * n + j] = dRR[2 * j + 0];
			x1[i] = 0.0;
		}
		// pick a random value for next x0
		CppAD::uniform_01(n, x0);
	}
}
void time_partial_hessian(size_t size, size_t repeat)
{	size_t J = size;
	size_t full_step = false;
	VECTOR(double) hess(2 * J * 2 * J), x(0);
	repeat_newton_hessian(repeat, J, full_step, x, hess);
}
void time_full_hessian(size_t size, size_t repeat)
{	size_t J = size;
	size_t full_step = true;
	VECTOR(double) hess(2 * J * 2 * J), x(0);
	repeat_newton_hessian(repeat, J, full_step, x, hess);
}
// END_REPEAT_NEWTON_HESSIAN
void center_title(size_t line_width, const std::string& title)
{	size_t skip = (line_width - title.size()) / 2;
	std::cout << "\n" << std::string(skip, ' ') << title << "\n";
	return;
}
// ---------------------------------------------------------------------------
int main(void)
{	bool ok       = true;
	double eps99  = 99.0 * std::numeric_limits<double>::epsilon();
	size_t J      = 10;
	size_t n      = 2 * J;
	size_t repeat = 1;
	bool   full_step, reverse;
	VECTOR(double) x(n);
	// ----------------------------------------------------------------------
	// check gradient results for a random x
	VECTOR(double) grad_kedem(n), grad_forward(n), grad_reverse(n);
	CppAD::uniform_01(n, x);
	repeat_kedem_gradient(repeat, J, x, grad_kedem);
	reverse = false;
	repeat_newton_gradient(repeat, J, reverse, x, grad_forward);
	reverse = true;
	repeat_newton_gradient(repeat, J, reverse, x, grad_reverse);
	for(size_t i = 0; i < n; i++)
	{	ok &= CppAD::NearEqual(grad_kedem[i], grad_forward[i], eps99, eps99);
		ok &= CppAD::NearEqual(grad_kedem[i], grad_reverse[i], eps99, eps99);
	}
	if( ok )
		std::cout << "Gradient values: OK\n";
	else
		std::cout << "Gradient values: Error\n";
	// -----------------------------------------------------------------------
	// check Hessian results for a random x
	VECTOR(double) hess_kedem(n*n), hess_partial(n*n), hess_full(n*n);
	CppAD::uniform_01(n, x);
	repeat_kedem_hessian(repeat, J, x, hess_kedem);
	full_step = false;
	repeat_newton_hessian(repeat, J, full_step, x, hess_partial);
	full_step = true;
	repeat_newton_hessian(repeat, J, full_step, x, hess_full);
	for(size_t i = 0; i < n*n; i++)
	{	ok &= CppAD::NearEqual(hess_kedem[i], hess_partial[i], eps99, eps99);
		ok &= CppAD::NearEqual(hess_kedem[i], hess_full[i],    eps99, eps99);
	}
	if( ok )
		std::cout << "Hessian values: OK\n";
	else
		std::cout << "Hessian values: Error\n";
	// -----------------------------------------------------------------------
	double time_min         = 1.0;
	size_t num_size         = 4;
	size_t size_width       = 5;
	size_t time_width       = 15;
	size_t time_precision   = 2;
	size_t ratio_width      = 10;
	size_t ratio_precision  = 5;
	VECTOR(std::string) name(3);
	using std::cout;
	using std::setw;
	using std::setprecision;
	cout << std::fixed;
	// -----------------------------------------------------------------------
	// Gradient Times
	VECTOR(size_t) grad_size(num_size);
	MATRIX(double) grad_time(num_size, 3);
	// title
	std::string  title = "Gradient Times in Milliseconds";
	size_t line_width  = size_width + 3 * time_width;
	center_title(line_width, title);
	// headings
	name[0] = "Kedem";
	name[1] = "Forward";
	name[2] = "Reverse";
	cout << setw(size_width) << "n";
	for(size_t j = 0; j < 3; j++)
		cout << setw(time_width) << name[j];
	cout << "\n";
	for(size_t i = 0; i < num_size; i++)
	{	J = 10 * (i + 1);
		grad_size[i]   = J;
		grad_time(i,0) = CppAD::time_test(time_kedem_gradient,   time_min, J);
		grad_time(i,1) = CppAD::time_test(time_forward_gradient, time_min, J);
		grad_time(i,2) = CppAD::time_test(time_reverse_gradient, time_min, J);
		cout << setw(size_width) << 2 * J;
		for(size_t j = 0; j < 3; j++)
		{	cout << setw(time_width) << setprecision(time_precision);
			cout << 1e3 * grad_time(i, j);
		}
		cout << "\n";
	}
	// -----------------------------------------------------------------------
	// Gradient Ratio
	VECTOR(std::string) factor_str(3);
	factor_str[0] = "n*n";
	factor_str[1] = "n*n";
	factor_str[2] = "n";
	MATRIX(double) factor(num_size, 3);
	for(size_t i = 0; i < num_size; i++)
	{	n = 2 * grad_size[i];
		factor(i,0) = double(n * n);
		factor(i,1) = double(n * n);
		factor(i,2) = double(n);
	}
	// title
	title      = "Gradient Milliseconds Divided by Factor";
	line_width = 10 + 10 + num_size * ratio_width;
	center_title(line_width, title);
	// headings
	cout << setw(10) << "Factor";
	cout << setw(10) << "n";
	for(size_t i = 0; i < num_size; i++)
		cout << setw(ratio_width) << 2 * grad_size[i];
	cout << "\n";
	// table values
	for(size_t j = 0; j < 3; j++)
	{	cout << setw(10) << factor_str[j];
		cout << setw(10) << name[j];
		for(size_t i = 0; i < num_size; i++)
		{	cout << setw(ratio_width) << setprecision(ratio_precision);
			cout << 1e3 * grad_time(i,j) / factor(i,j);
		}
		cout << "\n";
	}
	// -----------------------------------------------------------------------
	// Hessian Times
	VECTOR(size_t) hess_size(num_size);
	MATRIX(double) hess_time(num_size, 3);
	// title
	title       = "Hessian Times in Milliseconds";
	line_width  = size_width + 3 * time_width;
	center_title(line_width, title);
	// headings
	name[0] = "Kedem";
	name[1] = "Partial";
	name[2] = "Full";
	cout << setw(size_width) << "n";
	for(size_t j = 0; j < 3; j++)
		cout << setw(time_width) << name[j];
	cout << "\n";
	// table values
	for(size_t i = 0; i < num_size; i++)
	{	J = 10 * (i + 1);
		hess_size[i]   = J;
		hess_time(i,0) = CppAD::time_test(time_kedem_hessian,  time_min, J);
		hess_time(i,1) = CppAD::time_test(time_partial_hessian, time_min, J);
		hess_time(i,2) = CppAD::time_test(time_full_hessian,    time_min, J);
		cout << setw(size_width) << 2 * J;
		for(size_t j = 0; j < 3; j++)
		{	cout << setw(time_width) << setprecision(time_precision);
			cout << 1e3 * hess_time(i, j);
		}
		cout << "\n";
	}
	// -----------------------------------------------------------------------
	// Hessian Ratio
	factor_str[0] = "n*n*(n+3)/2";
	factor_str[1] = "n*(n+1)";
	factor_str[2] = "n*(n+1)";
	for(size_t i = 0; i < num_size; i++)
	{	n = 2 * hess_size[i];
		factor(i,0) = double(n * n * (n+3)) / 2.0;
		factor(i,1) = double(n * (n+1));
		factor(i,2) = double(n * (n+1));
	}
	// title
	title      = "Hessian Milliseconds Divided by Factor";
	line_width = 15 + 10 + num_size * ratio_width;
	center_title(line_width, title);
	// headings
	cout << setw(15) << "Factor";
	cout << setw(10) << "n";
	for(size_t i = 0; i < num_size; i++)
		cout << setw(ratio_width) << 2 * hess_size[i];
	cout << "\n";
	// table values
	for(size_t j = 0; j < 3; j++)
	{	cout << setw(15) << factor_str[j];
		cout << setw(10) << name[j];
		for(size_t i = 0; i < num_size; i++)
		{	cout << setw(ratio_width) << setprecision(ratio_precision);
			cout << 1e3 * hess_time(i,j) / factor(i,j);
		}
		cout << "\n";
	}
	// -----------------------------------------------------------------------
	if( ok )
		return 0;
	return 1;
}
