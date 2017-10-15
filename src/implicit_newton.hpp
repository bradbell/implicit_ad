# ifndef IMPLICIT_NEWTON_HPP
# define IMPLICIT_NEWTON_HPP

# include <cppad/example/cppad_eigen.hpp>
# include "utility.hpp"
/*
$begin implicit_newton$$
$spell
	num
	dw
	bool
	CppAD
	cppad
	rk
	xk
	const
	Taylor
$$

$section Newton Step Method for Derivatives of Implicit Functions$$

$head Syntax$$
$codei%implicit_newton %newton_ad%(
	%full_step%, %num_step%, %aL_fun%, %F_fun%, %solve%
)
%$$
$icode%rk% = %newton_ad%.Forward(%k%, %xk%)
%$$
$icode%dw% = %newton_ad%.Reverse(%q%, %w%)
%$$

$head Purpose$$
Given a function
$latex L : \B{R}^{n \times m} \rightarrow \B{R}^m$$,
we define the implicit function
$latex Y : \B{R}^n \rightarrow \B{R}^m$$ by
$latex L[ x , Y(x) ] = 0$$.
The partial $latex L_y [ x , Y(x) ]$$
is assumed to be invertible for all $latex x$$.
We define the reduced function
$latex R(x) = F[ x , Y(x) ]$$ where
$latex F : \B{R}^{n \times m} \rightarrow \B{R}^p$$.
The object $icode newton_ad$$ can be used to compute
derivatives of the reduced function $latex R(x)$$.
Note that in the special case where $latex F(x, y) = y$$,
$latex R(x) = Y(x)$$.

$head full_step$$
This argument has prototype
$codei%
	bool %full_step%
%$$
If it is true, full Newton steps are used.
If it is false, partial Newton steps are used.

$head num_step$$
This argument has prototype
$codei%
	size_t %num_step%
%$$
It is the number of Newton steps in the method.

$head aL_fun$$
This argument has prototype
$codei%
	CppAD::ADFun< CppAD::AD<double> >& %aL_fun%
%$$
and is the CppAD function object corresponding to $latex L(x, y)$$.
Note that a $code CppAD::ADFun<double>$$ object could be used
for the partial Newton step method.
Note that this function is not $code const$$ because
its forward mode is used for calculations.
Upon return,
there are no forward mode coefficient left in $icode aL_fun$$.

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
	const VECTOR(%Scalar%)&  %x%
	const VECTOR(%Scalar%)&  %y%
%$$
The return value has prototype
$codei%
	VECTOR(%Scalar%) %b%
%$$
This returns the value of $latex L_y (x, y)$$ for subsequent
calls to $icode%solve%.linear%$$.
Only the elements of $latex L_y (x, y)$$ that depend on $latex (x, y)$$
need be included in the vector $latex b$$.
The type $icode Scalar$$ is either $code double$$ or $code CppAD::AD<double>$$.

$subhead solve.linear$$
In the syntax
$codei%
	%u% = %solve%.linear(%b%, %v%)
%$$
the arguments $icode b$$, $icode v$$ and the return value
$icode u$$ have prototypes
$codei%
	const VECTOR(CppAD::AD<double>)& %b%
	const VECTOR(CppAD::AD<double>)& %v%
	const VECTOR(CppAD::AD<double>)  %u%
%$$
where both vectors have size $icode m$$.
The argument $icode b$$ has prototype
$codei%
	VECTOR(CppAD::AD<double>)& %b%
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
the previous call to $code forward$$ must have for
order $icode%k% - 1%$$ or greater.

$head rk$$
The return value $icode rk$$ has prototype
$codei%
	VECTOR(double) %rk%
%$$
its size is $icode m$$ and it is the $th k$$ order Taylor coefficient
for $icode R(x)$$.

$head q$$
This argument has prototype
$codei%
	size_t %q%
%$$
and is the number of Taylor coefficients in the function
that we are differentiating.

$head w$$
This vector has prototype
$codei%
	const VECTOR(double)& %w%
%$$
and its size is $icode%p%*%q%$$.
For $latex i = 0 , \ldots , p-1$$ and $latex k = 0 , \ldots , q-1$$,
let $latex r_i^k$$ denote the
$th k$$ Taylor coefficients for $latex R_i (x)$$
The function we are differentiating is
$latex \[
	W(r) = \sum_{i=0}^{p-1} \sum_{k=0}^{q-1} w_{i * q + k} r_i^k
\] $$

$head dw$$
This return value has prototype
$codei%
	VECTOR(double) %dw%
%$$
and its size is $icode%n%*%q%$$.
For $latex i = 0 , \ldots , n-1$$ and $latex k = 0 , \ldots , q-1$$,
let $latex X_i^k$$ denote the $th k$$ Taylor coefficients for $latex x$$.
The value $icode%dw%[%i%*%q%+%k%]%$$ is the partial of
$latex W[r(X)]%$$ with respect to $latex X_i^k$$.

$children%
	src/implicit_newton.cpp
%$$
$head Example$$
The routine $cref test_circle_implicit_newton$$
is a simple example and test using this class.
The routine $cref test_control_reduced_objective$$ is a control problem
example and test using this class.

$end
*/
// ----------------------------------------------------------------------------
// class definition
// ----------------------------------------------------------------------------
template <class Solve>
class implicit_newton {
private:
    // Newton step approximation for Y(x)`
	CppAD::ADFun<double> Y_fun_;
    // corresponds to F(x, y)
	CppAD::ADFun<double> F_fun_;
	// solves for y in L(x, y) = 0
	Solve                solve_;
	// number of orders currently stored in Y_fun_ and F_fun_.
	size_t               n_order_;
	// number of elements used to specify a value for L_y
	size_t               nb_;
public:
	implicit_newton(
		bool                               full_step  ,
		size_t                             num_step   ,
		CppAD::ADFun< CppAD::AD<double> >& aL_fun_    ,
		const CppAD::ADFun<double>&        F_fun_     ,
		const Solve&                       solve
	);
	VECTOR(double) Forward(size_t k, const VECTOR(double)& xk);
	VECTOR(double) Reverse(size_t k, const VECTOR(double)& w);
};
// ----------------------------------------------------------------------------
// class implementation
// ----------------------------------------------------------------------------
template <class Solve>
implicit_newton<Solve>::implicit_newton(
	bool                                    full_step  ,
	size_t                                  num_step   ,
	CppAD::ADFun< CppAD::AD<double>>&       aL_fun     ,
	const CppAD::ADFun<double>&             F_fun      ,
	const Solve&                            solve      )
{	typedef CppAD::AD<double> adouble;
	//
	assert( aL_fun.Domain() > aL_fun.Range() );
	assert( F_fun.Domain() == aL_fun.Domain() );
	size_t m = aL_fun.Range();
	size_t n = aL_fun.Domain() - m;
	//
	// F_fun_, solve_, n_order_
	F_fun_   = F_fun;
	solve_   = solve;
	n_order_ = 0;
	// ---------------------------------------------------------------------
	// Y_fun_
	// ---------------------------------------------------------------------
	// x
	VECTOR(double) x(n);
	for(size_t i = 0; i < n; i++)
		x[i] = 0.0;
	//
	// y^0
	VECTOR(double) y = solve_.function(x);
	assert( size_t( y.size() ) == m );
	//
	// L_y ( x, y^0 )
	VECTOR(double) b = solve_.derivative(x, y);
	//
	// set nb_
	nb_ = size_t( b.size() );
	//
	// start recording for Y_fun_(x, y, b)
	VECTOR(adouble) axyb0(n + m + nb_);
	for(size_t i = 0; i < n; i++)
		axyb0[i] = x[i];
	for(size_t i = 0; i < m; i++)
		axyb0[n + i] = y[i];
	for(size_t i = 0; i < nb_; i++)
		axyb0[n + m + i] = b[i];
	CppAD::Independent(axyb0);
	//
	// num_step Newton steps
	VECTOR(adouble) ax(n), ay(m), axy(n + m), axyb(n + m + nb_);
	axyb = axyb0;
	for(size_t j = 0; j < num_step; j++)
	{	//
		// L (x, y^j )
		for(size_t i = 0; i < n + m; i++)
			axy[i] = axyb[i];
		VECTOR(adouble) aL = aL_fun.Forward(0, axy);
		//
		// full:    L_y ( x, y^j ) step = - L ( x, y^j )
		// partial: L_y ( x, y^0 ) step = - L ( x, y^j )
		VECTOR(adouble) ab( nb_ );
		for(size_t i = 0; i < nb_; i++)
			ab[i] = axyb[n + m + i];
		VECTOR(adouble) astep = - solve_.linear(ab, aL);
		//
		// full:    y^{j+1} = y^j - L_y ( x , y^j )^{-1} L ( x, y^j )
		// partial: y^{j+1} = y^j - L_y ( x , y^0 )^{-1} L ( x, y^j )
		for(size_t i = 0; i < m; i++)
			axyb[n + i] = axyb[n + i] + astep[i];
		//
		if( full_step && (j + 1 < num_step) )
		{	// L_y (x, y^{j+1} )
			for(size_t i = 0; i < n; i++)
				ax[i] = axyb[i];
			for(size_t i = 0; i < m; i++)
				ay[i] = axyb[n + i];
			ab = solve_.derivative(ax, ay);
			for(size_t i = 0; i < nb_; i++)
				axyb[n + m + i] = ab[i];
		}
	}
	//
	// end recording of Y_fun_
	for(size_t i = 0; i < m; i++)
		ay[i] = axyb[n + i];
	Y_fun_.Dependent(axyb0, ay);
	//
	// clear the Taylor coefficients stored in aL_fun
	aL_fun.capacity_order(0);
}
template <class Solve> VECTOR(double)
implicit_newton<Solve>::Forward(size_t k, const VECTOR(double)& xk)
{	size_t m = Y_fun_.Range();
	size_t n = F_fun_.Domain() - m;
	size_t p = F_fun_.Range();
	//
	VECTOR(double) yk(m), xyk(n + m), bk(nb_), xybk(n + m + nb_), rk(p);
	if( k == 0 )
	{	// compute y0
		yk = solve_.function(xk);
		//
		bk = solve_.derivative(xk, yk);
		//
		// store zero order forward results in Y_fun_
		join_vector(xyk,  xk, yk);
		join_vector(xybk, xyk, bk);
		Y_fun_.Forward(k, xybk);
		//
		// get and store zero order forward results in F_fun_
		rk = F_fun_.Forward(k, xyk);
		//
		// number of orders currently stored
		n_order_ = 1;
		//
		// return result of zero order forward for F_fun_
		return rk;
	}
	assert( n_order_ >= k );
	for(size_t i = 0; i < m; i++)
		yk[i] = 0.0;
	for(size_t i = 0; i < nb_; i++)
		bk[i] = 0.0;
	//
	// compute and store k-th order result for Y_fun_
	join_vector(xyk, xk, yk);
	join_vector(xybk, xyk, bk);
	yk = Y_fun_.Forward(k, xybk);
	//
	// compute and store k-th order result for F_fun_
	join_vector(xyk, xk, yk);
	rk = F_fun_.Forward(k, xyk);
	//
	// number of orders currently stored
	n_order_ = k + 1;
	//
	return rk;
}
template <class Solve> VECTOR(double)
implicit_newton<Solve>::Reverse(size_t q, const VECTOR(double)& w)
{	// check that we have computed q-1 order forward mode
	assert(n_order_ >= q);
	//
	size_t m = Y_fun_.Range();
	size_t n = F_fun_.Domain() - m;
	//
	// check size of w
	assert( size_t(w.size()) == q * F_fun_.Range() );
	//
	// reverse F_fun
	VECTOR(double) dw_dxy = F_fun_.Reverse(q, w);
	assert( size_t( dw_dxy.size() ) == (n + m) * q );
	//
	// F_y [ x, Y(x) ] * Y'(x)
	VECTOR(double) dw_dy(m * q);
	for(size_t i = 0; i < m; i++)
	{	for(size_t k = 0; k < q; k++)
			dw_dy[i * q + k] = dw_dxy[ (n + i) * q + k];
	}
	VECTOR(double) dw_dxyb = Y_fun_.Reverse(q, dw_dy);
	assert( size_t( dw_dxyb.size() ) == (n + m + nb_) * q );
	//
	// F_x [x , Y(x)] + F_y [ x , Y(x) ] * Y'(x)
	VECTOR(double) dw_dx(n * q);
	for(size_t i = 0; i < n; i++)
	{	for(size_t k = 0; k < q; k++)
			dw_dx[i * q + k] = dw_dxy[i * q + k] + dw_dxyb[i * q + k ];
	}
	//
	return dw_dx;
}
# endif
