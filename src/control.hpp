# ifndef CONTROL_HPP
# define CONTROL_HPP

# include "utility.hpp"

namespace control { // BEGIN_CONTROL_NAMESPACE
/*
$begin vector_matrix$$
$spell
	xy
	vec
$$

$section Conversions Between Control Vectors and Matrices$$

$head xy_vec2mat$$
$srccode%cpp% */
template <class Scalar> void xy_vec2mat(
const VECTOR(Scalar)& xy_vec , MATRIX(Scalar)& x_mat, MATRIX(Scalar)& y_mat )
{	size_t J = size_t( x_mat.cols() );
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 2; i++)
			x_mat(i, j) = xy_vec[ j * 2 + i ];
		for(size_t i = 0; i < 4; i++)
			y_mat(i, j) = xy_vec[ 2 * J + j * 4 + i ];
	}
}
/* %$$
$end
*/


/*
-----------------------------------------------------------------------------
$begin objective$$

$section Computes the Control Objective Function$$

$head Syntax$$
$icode%F% = control::objective(%delta_t%, %q%, %x%, %y%)%$$

$head Prototype$$
$srcfile%control.hpp
	%0%// BEGIN_OBJECTIVE_PROTOTYPE%// END_OBJECTIVE_PROTOTYPE%
1%$$

$head Definition$$
The objective
$latex F : \B{R}^{2 \times J} \times \B{R}^{4 \times J} \rightarrow \B{R} $$
is defined by
$latex \[
F (x, y)
=
\frac{1}{2} \sum_{i=0}^3 q_i y_{i,J-1}^2
+
\frac{\Delta t}{4} \left(
	x_{0,0}^2 + x_{1,0}^2 + x_{0,J-1}^2 + x_{1,J-1}^2
	+
	2 \sum_{j=1}^{J-2} \left( x_{0,j}^2 + x_{1,j}^2 \right)
\right)
\] $$

$head Example$$
$srcfile%control.cpp%
	0%// BEGIN_TEST_OBJECTIVE%// END_TEST_OBJECTIVE%$$

$end
*/
// BEGIN_OBJECTIVE_PROTOTYPE
template <class Scalar> Scalar
objective(
	const Scalar&         delta_t ,
	const VECTOR(Scalar)& q       ,
	const MATRIX(Scalar)& x       ,
	const MATRIX(Scalar)& y       )
// END_OBJECTIVE_PROTOTYPE
{	assert( q.size() == 4 );
	assert( x.rows() == 2 );
	assert( y.rows() == 4 );
	assert( x.cols() == y.cols() );
	size_t J = size_t( x.cols() );
	//
	Scalar sum = Scalar(0.0);
	sum += x(0,0) * x(0,0)     + x(1,0) * x(1,0);
	sum += x(0,J-1) * x(0,J-1) + x(1,J-1) * x(1,J-1);
	for(size_t j = 1; j < J-1; j++)
		sum += Scalar(2.0) * ( x(0,j) * x(0,j) + x(1,j) * x(1,j) );
	sum = delta_t * sum / Scalar(4.0);
	//
	for(size_t i = 0; i < 4; i++)
		sum += q[i] * y(i, J-1) * y(i, J-1) / Scalar(2.0);
	//
	return sum;
}
/*
-------------------------------------------------------------------------------
$begin rec_objective$$
$spell
	CppAD
$$

$section Record the Control Objective$$

$head Syntax$$
$codei%control::rec_objective(%F_fun%, %J%, %delta_t%, %q%)%$$

$head Prototype$$
$srcfile%control.hpp
	%0%// BEGIN_REC_OBJECTIVE_PROTOTYPE%// END_REC_OBJECTIVE_PROTOTYPE%
1%$$

$head F_fun$$
The $cref objective$$ function
$latex F : \B{R}^{2 \times J} \times \B{R}^{4 \times J} \rightarrow \B{R} $$
is stored in $icode F_fun$$.
We use $icode u$$ to denote an $icode F_fun$$
argument vector (size $icode%6%*%J%$$)
and $icode v$$ a result vector (size one); e.g.,
$codei%
	%v% = %F_fun%.Forward(0, %u%)
%$$

$subhead J$$
is the number of time points in the control problem.

$subhead delta_t$$
is the step size between time points.

$head q$$
Is the objective weights for the components of the final state vector
in the control problem.

$subhead x$$
For $latex i = 0, 1$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	u_{j * 2 + i} = x_{i,j}
\] $$

$subhead y$$
For $latex i = 0, \ldots, 3$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	u_{2 * J + j * 4 + i} = y_{i,j}
\] $$

$subhead L$$
For $latex i = 0, \ldots, 3$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	v_{j * 4 + i} = L_{i,j}
\] $$

$end
*/
// BEGIN_REC_OBJECTIVE_PROTOTYPE
template <class Base> void
rec_objective(
	CppAD::ADFun<Base>&             F_fun   ,
	size_t                          J       ,
	const Base                      delta_t ,
	const VECTOR(Base)&             q       )
// END_REC_OBJECTIVE_PROTOTYPE
{	typedef CppAD::AD<Base> Scalar;
	//
	Scalar adelta_t = Scalar(delta_t);
	//
	VECTOR(Scalar) aq(4);
	for(size_t i = 0; i < 4; i++)
		aq[i] = Scalar( q[i] );
	//
	// declare independent variable vector
	VECTOR(Scalar) axy_vec(6 * J);
	for(size_t j = 0; j < 6 * J; j++)
		axy_vec[j] = 0.0;
	CppAD::Independent(axy_vec);
	//
	// convert to matrices in column major order
	MATRIX(Scalar) ax(2, J), ay(4, J);
	xy_vec2mat(axy_vec, ax, ay);
	//
	// Evaluate F(x, y)
	VECTOR(Scalar) aF_vec(1);
	aF_vec[0] = objective(adelta_t, aq, ax, ay);
	//
	// store the recording in F_fun
	F_fun.Dependent(axy_vec, aF_vec);
	return;
}
/*
-------------------------------------------------------------------------------
$begin constraint$$
$spell
$$

$section Computes the Control Constraint Function$$

$head Syntax$$
$icode%L% = control::constraint(%delta_t%, %p%, %x%, %y%)%$$

$head Prototype$$
$srcfile%control.hpp
	%0%// BEGIN_CONSTRAINT_PROTOTYPE%// END_CONSTRAINT_PROTOTYPE%
1%$$

$head Definition$$
The constraint equation is $latex L(x, y) = 0$$ where
$latex
L : \B{R}^{2 \times J} \times \B{R}^{4 \times J} \rightarrow \B{R}^{4 \times J}
$$
is defined by
$latex \[
	L_{i,0} (x, y) =  y_{i,0} -  p_i
\] $$
for $latex i = 0, \ldots , 3$$.
$latex \[
\begin{array}{rcl}
	a_j (y)        & = & [ ( y_{0,j} + 1 )^2 + y_{1,j}^2 ]^{-3/2} - 1
	\\
	L_{0,j} (x, y) & = & y_{2,j}
	\\
	L_{1,j} (x, y) & = & y_{3,j}
	\\
	L_{2,j} (x, y) & = & 2 y_{3,j} - (1 + y_{0,j}) a_j(y) + x_{0,j}
	\\
	L_{3,j} (x, y) & = & - 2 y_{2,j} - y_{1,j} a_j(y) + x_{1,j}
\end{array}
\] $$
for $latex j = 1, \ldots , J - 1$$.

$head Example$$
$srcfile%control.cpp%
	0%// BEGIN_TEST_CONSTRAINT%// END_TEST_CONSTRAINT%$$

$end
*/
// BEGIN_CONSTRAINT_PROTOTYPE
template <class Scalar> MATRIX(Scalar)
constraint(
	const Scalar&         delta_t ,
	const VECTOR(Scalar)& p       ,
	const MATRIX(Scalar)& x       ,
	const MATRIX(Scalar)& y       )
// END_CONSTRAINT_PROTOTYPE
{	assert( p.size() == 4 );
	assert( x.rows() == 2 );
	assert( y.rows() == 4 );
	assert( x.cols() == y.cols() );
	size_t J    = size_t( x.cols() );
	Scalar one  = Scalar(1.0);
	Scalar two  = Scalar(2.0);
	//
	MATRIX(Scalar) L(4,J);
	for(size_t i = 0; i < 4; i++)
		L(i,0) = y(i, 0) - p[i];
	for(size_t j = 1; j < J; j++)
	{	Scalar x0 = x(0,j-1);
		Scalar x1 = x(1,j-1);
		Scalar y0 = y(0,j-1);
		Scalar y1 = y(1,j-1);
		Scalar y2 = y(2,j-1);
		Scalar y3 = y(3,j-1);
		//
		Scalar r  = sqrt( (y0 + one) * (y0 + one) + y1 * y1 );
		Scalar a  = one / (r * r * r) - one;
		//
		Scalar dydt_0 = y2;
		Scalar dydt_1 = y3;
		Scalar dydt_2 = + two * y3 - (one + y0) * a + x0;
		Scalar dydt_3 = - two * y2 -         y1 * a + x1;
		//
		L(0,j)    = y(0,j) - y0 - delta_t * dydt_0;
		L(1,j)    = y(1,j) - y1 - delta_t * dydt_1;
		L(2,j)    = y(2,j) - y2 - delta_t * dydt_2;
		L(3,j)    = y(3,j) - y3 - delta_t * dydt_3;
	}
	//
	return L;
}
/*
-------------------------------------------------------------------------------
$begin rec_constraint$$
$spell
	CppAD
$$

$section Record the Control Constraint as a CppAD Function Object$$

$head Syntax$$
$codei%control::rec_constraint(%L_fun%, %J%, %delta_t%, %p%)%$$

$head Prototype$$
$srcfile%control.hpp
	%0%// BEGIN_REC_CONSTRAINT_PROTOTYPE%// END_REC_CONSTRAINT_PROTOTYPE%
1%$$

$head L_fun$$
The $cref constraint$$ function
$latex
L : \B{R}^{2 \times J} \times \B{R}^{4 \times J} \rightarrow \B{R}^{4 \times J}
$$
is stored in $icode L_fun$$.
We use $icode u$$ to denote an $icode L_fun$$
argument vector (size $icode%6%*%J%$$)
and $icode v$$ a result vector (size $icode%4%*%J%$$); e.g.,
$codei%
	%v% = %L_fun%.Forward(0, %u%)
%$$

$subhead J$$
is the number of time points in the control problem.

$subhead delta_t$$
is the step size between time points.

$head p$$
Is the initial value for the state variable vector
in the control problem.

$subhead x$$
For $latex i = 0, 1$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	u_{j * 2 + i} = x_{i,j}
\] $$

$subhead y$$
For $latex i = 0, \ldots, 3$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	u_{2 * J + j * 4 + i} = y_{i,j}
\] $$

$subhead L$$
For $latex i = 0, \ldots, 3$$ and $latex j = 0, \ldots , J-1$$,
$latex \[
	v_{j * 4 + i} = L_{i,j}
\] $$

$head Example$$

$subhead ADFun<double>$$
The following example records the constraint as an
$code ADFun<double>$$ object.
$srcfile%control.cpp%
	0%// BEGIN_TEST_REC_CONSTRAINT%// END_TEST_REC_CONSTRAINT%$$

$subhead ADFun< AD<double> >$$
The following example records the constraint as an
$code ADFun< AD<double> >$$ object.
$srcfile%control.cpp%
	0%// BEGIN_TEST_AD_REC_CONSTRAINT%// END_TEST_AD_REC_CONSTRAINT%$$

$end
*/
// BEGIN_REC_CONSTRAINT_PROTOTYPE
template <class Base> void
rec_constraint(
	CppAD::ADFun<Base>&             L_fun   ,
	size_t                          J       ,
	const Base                      delta_t ,
	const VECTOR(Base)&             p       )
// END_REC_CONSTRAINT_PROTOTYPE
{	typedef CppAD::AD<Base> Scalar;
	//
	Scalar adelta_t = Scalar(delta_t);
	//
	VECTOR(Scalar) ap(4);
	for(size_t i = 0; i < 4; i++)
		ap[i] = Scalar( p[i] );
	//
	// declare independent variable vector
	VECTOR(Scalar) axy_vec(6 * J);
	for(size_t j = 0; j < 6 * J; j++)
		axy_vec[j] = 0.0;
	CppAD::Independent(axy_vec);
	//
	// convert to matrices in column major order
	MATRIX(Scalar) ax(2, J), ay(4, J);
	xy_vec2mat(axy_vec, ax, ay);
	//
	// Evaluate L(x, y)
	MATRIX(Scalar) aL = constraint(adelta_t, ap, ax, ay);
	//
	// convert to a vector in column major order
	VECTOR(Scalar) aL_vec(4 * J);
	for(size_t j = 0; j < J; j++)
	{	for(size_t i = 0; i < 4; i++)
			aL_vec[j * 4 + i] = aL(i, j);
	}
	//
	// store the recording in L_fun
	L_fun.Dependent(axy_vec, aL_vec);
	return;
}
/*
-------------------------------------------------------------------------------
$begin full_newton$$
$spell
	num_itr
	xy
	Taylor
	CppAD
	rcv
	nr
$$

$section Execute Full Newton Steps For Control Constraint$$

$head Syntax$$
$icode%num_itr% = control::full_newton(
%xy_out%, %xy_in%, %L_fun%, %criteria%, %max_itr%, %L_y%, %work%
)%$$

$head Prototype$$
$srcfile%control.hpp
	%0%// BEGIN_FULL_NEWTON_PROTOTYPE%// END_FULL_NEWTON_PROTOTYPE%
1%$$

$head Purpose$$
Let $latex ( x^0 , y^0 )$$ denote the value of $latex (x, y)$$
in $icode xy_in$$.
Let $latex y^k$$ denote the value of $latex y$$
after $icode k$$ Newton steps (not its $th k$$ order Taylor coefficient).
The $th k$$ Newton step is
$latex \[
	y^k = y^{k-1} - L_y ( x^{k-1} , y^{k-1} )^{-1} L ( x^{k-1}, y^{k-1} )
\] $$

$head xy_in$$
This vector has size $icode%n% + %m%$$.
Its first $icode n$$ components specify $icode x$$.
The other components specify the initial value for $icode y$$.
It is the value of $latex ( x^0 , y^0 )$$ for the first Newton step.

$head xy_out$$
This vector has size $icode%n% + %m%$$.
Its first $icode n$$ components specify $icode x$$.
The other components specify the initial value for $icode y$$.
It is the value of $latex ( x^k , y^k )$$ for the last Newton step.

$head L_fun$$
The operation sequence for the $cref constraint$$ function
$latex L : \B{R}^{2 J} \times \B{R}^{4 J} \rightarrow \B{R}^{4 J}$$
is stored in $icode L_fun$$.

$head criteria$$
This is the convergence criteria in terms of the Euclidean norm squared.
Convergence is accepted when $latex | L(x, y) |^2$$ is less than
$icode criteria$$.

$head max_itr$$
This is the maximum number of Newton steps to execute
(which must be greater than zero).
If $icode%criteria% = 0%$$,
$icode xy_out$$ will correspond to exactly $icode max_itr$$ Newton steps.

$head L_y$$
This argument is the value of $latex L_y (x, y)$$ at the
value of $latex (x, y)$$ corresponding to $icode xy_in$$.

$head work$$
This is a work vector used to reduce the work.
It must correspond to $icode L_y$$ and not be empty;
see $cref jac_constraint$$.

$head num_itr$$
Is the number of Newton steps executed.
If $icode%num_itr% == %max_itr%$$,
the convergence criteria may not be satisfied.

$head Example$$
$srcfile%control.cpp%
	0%// BEGIN_TEST_FULL_NEWTON%// END_TEST_FULL_NEWTON%$$

$end
*/
// BEGIN_FULL_NEWTON_PROTOTYPE
template <class Scalar> size_t
full_newton(
	VECTOR(Scalar)&                                    xy_out        ,
	const VECTOR(Scalar)&                              xy_in         ,
	CppAD::ADFun<Scalar>&                              L_fun         ,
	Scalar                                             criteria      ,
	size_t                                             max_itr       ,
	const CPPAD_SPARSE(Scalar)&                        L_y           ,
	CppAD::sparse_jac_work&                            work          )
// END_FULL_NEWTON_PROTOTYPE
{	assert( size_t( xy_in.size() ) == L_fun.Domain() );
	assert( 0 < max_itr );
	size_t m = L_fun.Range();
	size_t n = L_fun.Domain() - m;
	assert( L_y.nr() != 0 );
	//
	// initilize Newton iterate
	xy_out = xy_in;
	//
	// L(x, y)
	VECTOR(Scalar) L = L_fun.Forward(0, xy_out);
	//
	// |L(x, y)|^2
	Scalar norm_sq_L = norm_squared(L);
	//
	if( norm_sq_L < criteria )
		return 0;
	//
	CPPAD_SPARSE(Scalar) L_y_itr = L_y;
	for(size_t itr = 0; itr < max_itr; itr++)
	{
		// matrix is lower triangular, solve for Newton step
		VECTOR(Scalar) step(m);
		solve_lower_cppad(L_y_itr, step, L);
		//
		// modify y components of xy_out
		for(size_t i = 0; i < m; i++)
			xy_out[n + i] -= step[i];
		//
		// L(x, y)
		L = L_fun.Forward(0, xy_out);
		//
		// |L(x, y)|^2
		norm_sq_L = norm_squared(L);
		//
		bool last = itr + 1 == max_itr || norm_sq_L < criteria;
		if( ! last )
		{	// evaluate L_y (x, y)
			jac_constraint(L_y_itr, L_fun, xy_out, work);
		}
		//
		if( norm_sq_L < criteria )
			return itr + 1;
	}
	return max_itr;
}
/*
-----------------------------------------------------------------------------
$begin implicit_solver$$
$spell
	CppAD
	cppad
	const
	Kedem
$$

$section Control Problem Solver for Implicit Kedem or Newton Object$$

$head Syntax$$
$codei%control::implicit_solver %solve%(%L_fun%, %aL_fun%, %a%criteria%)
%$$
$icode%y% = %solve%.function(%x%)
%$$
$icode%b% = %solve%.derivative(%x%, %y%)
%$$
$icode%u% = %solve%.linear(%v%)%$$

$head Purpose$$
The object $icode solve$$ can be used with
$cref/implicit_kedem/implicit_kedem/solve/$$ to compute
derivatives of the function
$latex Y(x)$$ defined by $latex L[x, Y(x)] = 0$$
for the control problem.

$head L_fun$$
This argument has prototype
$codei%
	const CppAD::ADFun<double>& %L_fun%
%$$
and is the CppAD function object corresponding to the
$cref constraint$$ function.

$head aL_fun$$
This argument has prototype
$codei%
	const CppAD::ADFun< CppAD::AD<double> >& %aL_fun%
%$$
and is the CppAD function object corresponding to the
$cref constraint$$ function.
This function object can be empty, if
$cref/solve.derivative/implicit_solver/solve.derivative/$$
is note used with $icode Scalar$$ equal to $code CppAD::AD<double>$$.
An empty $icode aL_fun$$ object can be created with
$codei%
	CppAD::ADFun< CppAD::AD<double> > %aL_fun%;
%$$

$head criteria$$
This is the convergence criteria for $latex L[x, Y(x)] = 0$$.
To be specific, a value $latex y$$ is accepted if
the Euclidean norm squared
$latex | L(x, y) |^2$$ is less than $icode criteria$$.

This argument has prototype
$codei%
	const %Solve%& %solve%
%$$
The type $icode Solver$$ must support the default constructor
and the assignment operator.
It must also support the following operations:

$head solve.function$$
The argument $icode x$$ and the return value $icode y$$ have prototypes
$codei%
	const VECTOR(double)&  %x%
	VECTOR(double)         %y%
%$$
where $icode%x%.size() == %n%$$ and $icode%y%.size() == %m%$$.
The return value satisfies the relation $latex L(x, y) = 0$$.

$head solve.derivative$$
The arguments have prototypes
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
The type $icode Scalar$$ is either
$code double$$ or $code CppAD::AD<double>$$.

$head solve.linear$$
The arguments $icode b$$, $icode v$$ and the return value
$icode u$$ have prototypes
$codei%
	const VECTOR(%Scalar%)& %b%
	const VECTOR(%Scalar%)& %v%
	VECTOR(%Scalar%)        %u%
%$$
where both vectors have size $icode m$$.
The return value satisfies the equation
$latex \[
	u = L_y (x, y)^{-1} v
\] $$
where $latex L_y (x, y)$$ corresponds to $icode b$$.
The type $icode Scalar$$ is either
$code double$$ or $code CppAD::AD<double>$$.

$head Example$$
$srcfile%control.cpp%
0%// BEGIN_TEST_CONTROL_IMPLICIT_SOLVER%// END_TEST_CONTROL_IMPLICIT_SOLVER%$$

$end
*/
class implicit_solver {
private:
	// function object representing L(x, y)
	CppAD::ADFun<double>               L_fun_;
	CppAD::ADFun< CppAD::AD<double> > aL_fun_;
	//
	// convergence criteria
	double criteria_;
	//
	// sparisty pattern for Jacobian of L(x, y) w.r.t. y
	CppAD::sparse_rc< VECTOR(size_t) > L_y_pattern_;
	//
	// work space used by for_sparse_jac
	CppAD::sparse_jac_work work_;
	//
public:
	// default constructor necessary for implicit_kedem
	implicit_solver(void)
	{ }
	implicit_solver(
		CppAD::ADFun<double>&              L_fun     ,
		CppAD::ADFun< CppAD::AD<double> >& aL_fun    ,
		double                             criteria  )
	{	assert( L_fun.Domain() > L_fun.Range() );
		size_t m = L_fun.Range();
		size_t n = L_fun.Domain() - m;
		assert( aL_fun.size_var() == 0 || (aL_fun.Domain() == n + m) );
		assert( aL_fun.size_var() == 0 || (aL_fun.Range() == m) );
		//
		// L_fun_
		L_fun_ = L_fun;
		//
		// aL_fun_
		aL_fun_ = aL_fun;
		//
		// criteria_
		criteria_ = criteria;
		//
		// work_
		CPPAD_SPARSE(double) L_y;
		VECTOR(double) xy(n + m);
		for(size_t i = 0; i < n + m; i++)
			xy[i] = 0.0;
		jac_constraint(L_y, L_fun_, xy, work_);
		//
		// L_y_pattern_
		size_t nr = L_y.nr();
		size_t nc = L_y.nc();
		size_t nnz = L_y.nnz();
		assert( nr == m );
		assert( nc == m );
		CppAD::sparse_rc< VECTOR(size_t) > pattern(nr, nc, nnz);
		for(size_t k = 0; k < nnz; k++)
			pattern.set(k, L_y.row()[k], L_y.col()[k]);
		L_y_pattern_ = pattern;
	}
	VECTOR(double) function(const VECTOR(double)& x)
	{	size_t m = L_fun_.Range();
		size_t n = L_fun_.Domain() - m;
		assert( size_t( x.size() ) == n );
		//
		size_t max_itr   = 10;
		VECTOR(double) xy_in(n + m), xy_out;
		for(size_t i = 0; i < n; i++)
			xy_in[i] = x[i];
		for(size_t j = 0; j < m; j++)
			xy_in[n + j] = 0.0;
		//
		// L_y (x, y)
		CPPAD_SPARSE(double) L_y( L_y_pattern_ );
		jac_constraint(L_y, L_fun_, xy_in, work_);
		//
		// set xy_out
# ifndef NDEBUG
		size_t num_itr   = full_newton(
			xy_out, xy_in, L_fun_, criteria_, max_itr, L_y, work_
		);
		// abort if convergence criteria not satisfied
		assert( num_itr < max_itr );
# else
		full_newton(xy_out, xy_in, L_fun_, criteria_, max_itr, L_y, work_);
# endif
		//
		// extract y subvector of xy_out s.t. L(x, y) = 0
		VECTOR(double) y(m);
		for(size_t j = 0; j < m; j++)
			y[j] = xy_out[n + j];
		//
		return y;
	}
	VECTOR(double) derivative(const VECTOR(double)& x, const VECTOR(double)& y)
	{	VECTOR(double) xy(x.size() + y.size());
		join_vector(xy, x, y);
		CPPAD_SPARSE(double) L_y( L_y_pattern_ );
		jac_constraint(L_y, L_fun_, xy, work_);
		return L_y.val();
	}
	VECTOR(CppAD::AD<double>) derivative(
		const VECTOR(CppAD::AD<double>)& ax ,
		const VECTOR(CppAD::AD<double>)& ay )
	{	VECTOR(CppAD::AD<double>) axy(ax.size() + ay.size());
		join_vector(axy, ax, ay);
		CPPAD_SPARSE(CppAD::AD<double>) aL_y( L_y_pattern_ );
		jac_constraint(aL_y, aL_fun_, axy, work_);
		return aL_y.val();
	}
	template <class Scalar>
	VECTOR(Scalar) linear(const VECTOR(Scalar)& b, const VECTOR(Scalar)& v)
	{	size_t nnz = L_y_pattern_.nnz();
		size_t nr  = L_y_pattern_.nr();
		//
		assert( L_y_pattern_.nc() == L_y_pattern_.nr() );
		assert( size_t( b.size() ) == nnz );
		assert( size_t( v.size() ) == nr );
		//
		// set values in L_y_ to b
		CPPAD_SPARSE(Scalar) L_y( L_y_pattern_ );
		for(size_t k = 0; k < nnz; k++)
			L_y.set(k, b[k]);
		//
		// u = L_y (x, y)^{-1} v
		VECTOR(Scalar) u( nr );
		solve_lower_cppad(L_y, u, v);
		//
		return u;
	}
};

} // END_CONTROL_NAMESPACE

# endif
