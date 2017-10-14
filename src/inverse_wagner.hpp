# ifndef INVERSE_WAGNER_HPP
# define INVERSE_WAGNER_HPP

# include <cppad/example/cppad_eigen.hpp>
# include <Eigen/Dense>

/*
$begin inverse_wagner$$
$spell
	cppad_inv
	xk
	zk
	CppAD
	typedef
	Eigen
	const
	Taylor
$$
$latex \newcommand{\B}[1]{{\bf #1}}$$

$section Wagner Method for Derivatives of Inverse Functions$$

$head Syntax$$
$codei%inverse_wagner %inv_fun%(%fun%, %solve%)
%$$
$icode%xk% = %inv_fun%.forward(%k%, %zk%)
%$$

$head Purpose$$
Given a CppAD function object corresponding to
$latex f : \B{R}^n \rightarrow \B{R}^n$$,
we create an object that implements forward mode for the inverse function.
The derivative $latex f^{(1)} (x)$$ must be invertible
for all $latex x$$ used by this routine.

$head vector$$
The type $code inverse_wagner::vector$$ is defined to be
$codep
	Eigen::Matrix<double, Eigen::Dynamic, 1>
$$

$head fun$$
This argument has prototype
$codei%
	const CppAD::ADFun<double>& %fun%
%$$
and is the CppAD function object corresponding to $latex z = f(x)$$.

$head solve$$
This argument has prototype
$codei%
	const %Solve%& %solve%
%$$
The type $icode Solver$$ must support the default constructor
and the assignment operator.

$subhead solve.function$$
The class $icode Solve$$ must support the syntax
$codei%
	%x% = %solve%.function(%z%)
%$$
The argument $icode z$$ and the return value $icode x$$ have prototypes
$codei%
	const inverse_wagner::vector&  %z%
	inverse_wagner::vector         %x%
%$$
where both vectors have size $icode n$$.
The return value satisfies the relation $latex z = f(x)$$.

$subhead solve.derivative$$
The class $icode Solve$$ must support the syntax
$codei%
	%u% = %solve%.derivative(%v%)
%$$
The argument $icode v$$ and the return value $icode u$$ have prototypes
$codei%
	const inverse_wagner::vector& %v%
	inverse_wagner::vector        %u%
%$$
where both vectors have size $icode n$$.
The return value satisfies the relation $latex f^{(1)} (z) u  = v$$
where $icode z$$ corresponds to the previous
$codei%
	x = %solve%.function(%z%)
%$$

$head k$$
This is the order for this forward mode calculation and has prototype
$codei%
	size_t %k%
%$$.

$subhead zk$$
The argument $icode zk$$ has prototypes
$codei%
	const inverse_wagner::vector& %zk%
%$$
its size is $icode n$$ and it is the $th k$$ order Taylor coefficient
for $icode z$$.
The results of the calculation are stored internally and used
during higher order forward mode calculations.
If $icode%k% > 0%$$,
there must have been a previous call to $code forward$$ for
order $icode%k% - 1%$$.

$head xk$$
The return value $icode xk$$ has prototype
$codei%
	inverse_wagner::vector %xk%
%$$
its size is $icode n$$ and it is the $th k$$ order Taylor coefficient
for $icode x$$.

$head Method$$
Fix $latex z : \B{R} \rightarrow \B{R}^n$$ and define
$latex x : \B{R} \rightarrow \B{R}^n$$ by
$latex z(t) = f[ x(t) ] $$.
The $th k$$ order Taylor coefficients are
$latex x^k = x^{(k)} (0) / k !$$ and
$latex z^k = z^{(k)} (0) / k !$$.


$subhead Zero Order$$
The zero order coefficient $latex x^0$$ is found using
$codei%
	%x0% = %solve%.function(%z0%)
%$$

$subhead Higher Orders$$
It follows from the definitions that
$latex \[
	z^{(1)} (t) = f^{(1)} ( x^0 ) x^{(1)} (t)
\] $$
Taking more derivatives, we see that, for $latex k \geq 1$$,
$latex \[
	z^k = f^{(1)} ( x^0 ) x^k + H_k ( x^0 , \cdots , x^{k-1} )
\] $$
We use $th k$$ order forward mode on $icode fun$$,
with $latex x^k = 0$$, to determine
$latex H_k ( x^0 , \cdots , x^{k-1} )$$.
We then solve the equation
$latex \[
	f^{(1)} ( x^0 ) x^k = z^k - H_k ( x^0 , \cdots , x^{k-1} )
\] $$
to determine $latex x^k$$.

$end
*/
// ----------------------------------------------------------------------------
// class definition
// ----------------------------------------------------------------------------
template <class Solve>
class inverse_wagner {
public:
	typedef typename Eigen::Matrix<double, Eigen::Dynamic, 1>  vector;
private:
    // corresponds to z = f(x)
	CppAD::ADFun<double> fun_;

	// inverse solver
	Solve                 solve_;

	// number of orders currently stored in fun_.
	size_t                n_order_;

	// zero vector (effectively const)
	vector                zero_;
public:
	inverse_wagner(
		const CppAD::ADFun<double>& fun_    ,
		const Solve&                solve
	);
	vector forward(size_t k, const vector& zk);
};
// ----------------------------------------------------------------------------
// class implementation
// ----------------------------------------------------------------------------
template <class Solve>
inverse_wagner<Solve>::inverse_wagner(
	const CppAD::ADFun<double>& fun   ,
	const Solve&                solve )
{	size_t n = fun.Domain();
	//
	assert( fun.Range() == n );
	fun_     = fun;
	solve_   = solve;
	n_order_ = 0;
	zero_.resize(n);
	for(size_t j = 0; j < n; j++)
		zero_[j] = 0.0;
}
template <class Solve> typename inverse_wagner<Solve>::vector
inverse_wagner<Solve>::forward(size_t k, const vector& zk)
{	vector xk( fun_.Domain() );
	if( k == 0 )
	{	// compute x0
		xk = solve_.function(zk);
		// store zero order forward results in fun_
		fun_.Forward(k, xk);
		// number of orders currently stored
		n_order_ = 1;
		//
		return xk;
	}
	assert( n_order_ >= k );
	// Hk = k-th order forward result with xk = 0.0
	vector Hk = fun_.Forward(k, zero_);
	// solve for x^k
	xk = solve_.derivative(zk - Hk);
	// store zero order forward result for this order
	fun_.Forward(k, xk);
	// number of orders currently stored
	n_order_ = k + 1;
	//
	return xk;
}
// ----------------------------------------------------------------------------
# endif
