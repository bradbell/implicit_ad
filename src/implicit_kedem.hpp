# ifndef IMPLICIT_KEDEM_HPP
# define IMPLICIT_KEDEM_HPP

# include <cppad/example/cppad_eigen.hpp>
# include "utility.hpp"
/*
$begin implicit_kedem$$
$spell
	CppAD
	cppad
	rk
	xk
	const
	Taylor
	Kedem
	bool
$$

$section Kedem Method for Derivatives of Implicit Functions$$

$head Syntax$$
$codei%implicit_kedem %kedem_ad%(%L_fun%, %F_fun%, %solve%)
%$$
$icode%rk% = %kedem_ad%.Forward(%k%, %xk%, %store%)
%$$

$head Purpose$$
Given a function
$latex L : \B{R}^{n \times m} \rightarrow \B{R}^m$$,
we define the implicit function
$latex Y : \B{R}^n \rightarrow \B{R}^m$$ by
$latex L[ x , Y(x) ] = 0$$.
The partial $latex L_y [ x , Y(x) ]$$
must be invertible for all $latex x$$ used by $icode kedem_ad$$.
We define the reduced function
$latex R(x) = F[ x , Y(x) ]$$ where
$latex F : \B{R}^{n \times m} \rightarrow \B{R}^p$$.
The object $icode kedem_ad$$ can be used to compute
derivatives of the reduced function $latex R(x)$$.
Note that in the special case where $latex F(x, y) = y$$,
$latex R(x) = Y(x)$$.

$head L_fun$$
This argument has prototype
$codei%
	const CppAD::ADFun<double>& %L_fun%
%$$
and is the CppAD function object corresponding to $latex L(x, y)$$.

$head F_fun$$
This argument has prototype
$codei%
	const CppAD::ADFun<double>& %F_fun%
%$$
and is the CppAD function object corresponding to $latex F(x, y)$$.

$head solve$$
This argument has prototype
$codei%
	const %Solve%& %solve%
%$$
The type $icode Solver$$ must support the default constructor
and the assignment operator.
It must also support the following operations:

$subhead solve.function$$
In the syntax
$codei%
	%y% = %solve%.function(%x%)
%$$
the argument $icode x$$ and the return value $icode y$$ have prototypes
$codei%
	const VECTOR(double)&  %x%
	VECTOR(double)         %y%
%$$
where $icode%x%.size() == %n%$$ and $icode%y%.size() == %m%$$.
The return value satisfies the relation $latex L(x, y) = 0$$.

$subhead solve.derivative$$
In the syntax
$codei%
	%b% = %solve%.derivative(%x%, %y%)
%$$
the arguments have prototypes
$codei%
	const VECTOR(double)&  %x%
	const VECTOR(double)&  %y%
%$$
The return value has prototype
$codei%
	VECTOR(double) %b%
%$$
This returns the value of $latex L_y (x, y)$$ for subsequent
calls to $icode%solve%.linear%$$.
Only the elements of $latex L_y (x, y)$$ that depend on $latex (x, y)$$
need be included in the vector $latex b$$.

$subhead solve.linear$$
In the syntax
$codei%
	%u% = %solve%.linear(%b%, %v%)
%$$
the arguments $icode b$$, $icode v$$ and the return value
$icode u$$ have prototypes
$codei%
	const VECTOR(double)& %b%
	const VECTOR(double)& %v%
	const VECTOR(double)  %u%
%$$
where both vectors have size $icode m$$.
The argument $icode b$$ has prototype
$codei%
	VECTOR(double)& %b%
%$$
The return value satisfies the equation
$latex \[
	u = L_y (x, y)^{-1} v
\] $$
where $latex L_y (x, y)$$ corresponds to $icode b$$.

$head k$$
This is the order for this forward mode calculation and has prototype
$codei%
	size_t %k%
%$$
Note that computing the $th k$$ order, uses the results of the
previous lower order calculations which are internally stored
in the function objects.

$head xk$$
The argument $icode xk$$ has prototypes
$codei%
	const VECTOR(double)& %xk%
%$$
its size is $icode m$$ and it is the $th k$$ order Taylor coefficient
for $icode x$$.
The results of the calculation are stored internally and used
during higher order forward mode calculations.
If $icode%k% > 0%$$,
there must have been a previous call to $code forward$$ for
order $icode%k% - 1%$$.

$head store$$
This argument has prototype
$codei%
	bool %store%
%$$
If it is true, the results of this forward mode calculation is store
so that a future forward mode calculation of order $icode%k%+1%$$ is possible.
This requires an extra pass over the operation sequence for
$cref/L_fun/implicit_kedem/L_fun/$$.
Hence $icode store$$ should be false when a higher order calculation
will not be needed.

$head rk$$
The return value $icode rk$$ has prototype
$codei%
	VECTOR(double) %rk%
%$$
its size is $icode m$$ and it is the $th k$$ order Taylor coefficient
for $icode R(x)$$.

$head Method$$
Fix $latex x : \B{R} \rightarrow \B{R}^n$$ and define
$latex y : \B{R} \rightarrow \B{R}^m$$ by
$latex x(t) = Y( x(t) )$$.
The $th k$$ order Taylor coefficients are
$latex x^k = x^{(k)} (0) / k !$$ and
$latex y^k = y^{(k)} (0) / k !$$.

$subhead Zero Order$$
The zero order coefficient $latex y^0$$ is found using
$codei%
	%y0% = %solve%.function(%x0%)
%$$

$subhead Higher Orders$$
It follows from the definitions that
$latex \[
	0 = L_x ( x^0 , y^0 ) x^{(1)} (t) + L_y ( x^0 , y^0 ) y^{(1)} (t)
\] $$
Taking more derivatives we see that, for $latex k \geq 1$$,
$latex \[
0 = H_k ( x^0 , \cdots , x^k , y^0 , \cdots , y^{k-1} ) + L_y ( x^0 , y^0 ) y^k
\] $$
We use $th k$$ order forward mode on $icode L_fun$$,
with $latex y^k = 0$$, to determine
$latex \[
	z^k = H_k ( x^0 , \cdots , x^k , y^0 , \cdots , y^{k-1} )
\] $$
We then solve the equation
$latex \[
	0 = z^k + L_y ( x^0 , y^0 ) y^k
\] $$
to determine the proper value for $latex y^k$$.

$children%
	implicit_kedem.cpp
%$$
$head Example$$
The routine $cref test_circle_implicit_kedem$$ is a simple example
and test using this class.
The routine $cref test_control_reduced_objective$$ is a control problem
example and test using this class.

$end
*/
// ----------------------------------------------------------------------------
// class definition
// ----------------------------------------------------------------------------
template <class Solve>
class implicit_kedem {
private:
    // corresponds to L(x, y)
	CppAD::ADFun<double> L_fun_;
	//
    // corresponds to F(x, y)
	CppAD::ADFun<double> F_fun_;
	//
	// implicit_kedem solver
	Solve                solve_;
	//
	// number of orders currently stored in L_fun_ and F_fun_.
	size_t               n_order_;
	//
	// value of L_y (x, y) corresponding to previous zero order forward
	// with store equal to true
	VECTOR(double)       b_;
public:
	implicit_kedem(
		const CppAD::ADFun<double>& L_fun_    ,
		const CppAD::ADFun<double>& F_fun_    ,
		const Solve&                solve
	);
	VECTOR(double) Forward(size_t k, const VECTOR(double)& xk, bool store);
};
// ----------------------------------------------------------------------------
// class implementation
// ----------------------------------------------------------------------------
template <class Solve>
implicit_kedem<Solve>::implicit_kedem(
	const CppAD::ADFun<double>& L_fun ,
	const CppAD::ADFun<double>& F_fun ,
	const Solve&                solve )
{	assert( L_fun.Domain() > L_fun.Range() );
	assert( F_fun.Domain() == L_fun.Domain() );
	L_fun_   = L_fun;
	F_fun_   = F_fun;
	solve_   = solve;
	n_order_ = 0;
}
template <class Solve> VECTOR(double)
implicit_kedem<Solve>::Forward(size_t k, const VECTOR(double)& xk, bool store)
{	size_t m = L_fun_.Range();
	size_t n = L_fun_.Domain() - m;
	size_t p = F_fun_.Range();
	//
	VECTOR(double) yk(m), xyk(n + m), rk(p);
	if( k == 0 )
	{	// compute y0
		yk = solve_.function(xk);
		//
		// xyk = (xk, yk)
		join_vector(xyk, xk, yk);
		//
		// get and store zero order forward results in F_fun_
		rk = F_fun_.Forward(k, xyk);
		//
		if( ! store )
		{	n_order_ = k;
			return rk;
		}
		//
		// compute and store L_y(x, y) corresponding to this (x, y)
		b_ = solve_.derivative(xk, yk);
		//
		// store zero order forward results in L_fun_
		L_fun_.Forward(k, xyk);
		//
		// number of orders currently stored
		n_order_ = k + 1;
		//
		// return result of zero order forward for F_fun_
		return rk;
	}
	assert( n_order_ >= k );
	// zk = k-th order forward result for L with yk = 0.0
	for(size_t i = 0; i < m; i++)
		yk[i] = 0.0;
	join_vector(xyk, xk, yk);
	VECTOR(double) minus_zk = - L_fun_.Forward(k, xyk);
	//
	// solve for y^k
	yk = solve_.linear(b_, minus_zk);
	//
	// xyk = (xk, yk)
	join_vector(xyk, xk, yk);
	//
	// get and store forward results for this order in F_fun_
	rk = F_fun_.Forward(k, xyk);
	//
	if( ! store )
	{	n_order_ = k;
		return rk;
	}
	// store forward result for this order in L_fun_
	L_fun_.Forward(k, xyk);
	//
	// number of orders currently stored
	n_order_ = k + 1;
	//
	return rk;
}
# endif
