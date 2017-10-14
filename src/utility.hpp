# ifndef UTILITY_HPP
# define UTILITY_HPP

# include <cppad/example/cppad_eigen.hpp>
# include <Eigen/Dense>
# include <Eigen/SparseCore>

/*
-------------------------------------------------------------------------------
$begin notation$$
$spell
	hpp
$$

$section Notation$$
The file $code utility.hpp$$ includes the following notation:


$head Dynamic$$
$srccode%cpp% */
using Eigen::Dynamic;
/* %$$

$head SparseMatrix$$
$srccode%cpp% */
using Eigen::SparseMatrix;
/* %$$

$head VECTOR$$
$srccode%cpp% */
# define VECTOR(Scalar)  Eigen::Matrix<Scalar, Dynamic, 1>
/* %$$

$head MATRIX$$
$srccode%cpp% */
# define MATRIX(Scalar)  Eigen::Matrix<Scalar, Dynamic, Dynamic>
/* %$$

$head CPPAD_SPARSE$$
$srccode%cpp% */
# define CPPAD_SPARSE(Scalar) \
	CppAD::sparse_rcv< VECTOR(size_t) , VECTOR(Scalar) >
/* %$$
$end
*/

/*
------------------------------------------------------------------------------
$begin norm_squared$$
$spell
	hpp
$$

$section Norm Squared of a Vector$$
The file $code utility.hpp$$ includes the following function:

$srccode%cpp% */
template <class Scalar> Scalar
norm_squared(const VECTOR(Scalar)& v)
{	return v.squaredNorm(); }
/*%$$
$end
*/

/*
------------------------------------------------------------------------------
$begin join_vector$$
$spell
	hpp
$$

$section Join Two Vectors$$
The file $code utility.hpp$$ includes the following function:

$srccode%cpp% */
template <class Scalar> void join_vector(
	VECTOR(Scalar)& xy, const VECTOR(Scalar)& x, const VECTOR(Scalar)& y
)
{	size_t n = size_t( x.size() );
	size_t m = size_t( y.size() );
	for(size_t i = 0; i < n; i++)
		xy[i] = x[i];
	for(size_t i = 0; i < m; i++)
		xy[n + i] = y[i];
}
/*%$$
$end
*/

/*
------------------------------------------------------------------------------
$begin sparse_cppad2eigen$$
$spell
	CppAD
	Eigen
	cppad
	cols
	hpp
$$

$section Convert A CppAD Sparse Matrix to an Eigen Sparse Matrix$$

$head Syntax$$
$codei%#include "utility.hpp"
%$$
$codei%sparse_cppad2eigen(%sparse_cppad%, %sparse_eigen)%$$.

$head Prototype$$
$srcfile%utility.hpp%
0%// BEGIN_SPARSE_CPPAD2EIGEN_PROTOTYPE%// END_SPARSE_CPPAD2EIGEN_PROTOTYPE%
1%$$

$head sparse_cppad$$
Is the CppAD sparse matrix that we are converting to an Eigen sparse matrix.

$head sparse_eigen$$
On input, if $icode%sparse_eigen%.rows()%$$ or $icode%sparse_eigen%.cols()%$$
is zero, a new sparsity patter is allocated in $icode sparse_eigen$$.
Otherwise, it is assumed that the sparsity pattern in
$icode sparse_eigen$$ corresponds to a previous call to
$code sparse_cppad2eigen$$
with the same sparsity pattern in $code cppad_sparse$$.

$head Example$$
$srcfile%utility.cpp%
	0%// BEGIN_TEST_SPARSE_CPPAD2EIGEN%// END_TEST_SPARSE_CPPAD2EIGEN%$$

$end
*/
// BEGIN_SPARSE_CPPAD2EIGEN_PROTOTYPE
template <class Scalar> void
	sparse_cppad2eigen(
	const CPPAD_SPARSE(Scalar)& sparse_cppad  ,
	SparseMatrix<Scalar>&       sparse_eigen  )
// END_SPARSE_CPPAD2EIGEN_PROTOTYPE
{
	// extract sparse_cppad information
	size_t nr  = sparse_cppad.nr();
	size_t nc  = sparse_cppad.nc();
	size_t nnz = sparse_cppad.nnz();
	const VECTOR(size_t)& row( sparse_cppad.row() );
	const VECTOR(size_t)& col( sparse_cppad.col() );
	const VECTOR(Scalar)& val( sparse_cppad.val() );
	//
	bool new_allocate = sparse_eigen.rows() == 0 || sparse_eigen.cols() == 0;
	if( new_allocate )
	{	//
		// count the number of entries in each column of sparse_cppad
		VECTOR(size_t) col_count(nc);
		for(size_t j = 0; j < nc; j++)
			col_count[j] = 0;
		for(size_t k = 0; k < nnz; k++)
		{	size_t c = col[k];
			++col_count[c];
		}
		//
		// size of matrix and allocate memory
		nr = sparse_cppad.nr();
		nc = sparse_cppad.nc();
		sparse_eigen.resize(nr, nc);
		sparse_eigen.reserve( col_count );
	}
	//
	// set the elements in the matrix
	for(size_t k = 0; k < nnz; k++)
	{	size_t r = row[k];
		size_t c = col[k];
		Scalar v = val[k];
		if( new_allocate )
			sparse_eigen.insert(r, c) = v;
		else
			sparse_eigen.coeffRef(r, c) = v;
	}
	return;
}
/*
-------------------------------------------------------------------------------
$begin solve_lower_cppad$$
$spell
	CppAD
	nr
	nc
	cppad
	hpp
$$

$section Solve a CppAD Sparse Lower Triangular System$$

$head Syntax$$
$codei%#include "utility.hpp"
%$$
$icode%msg% = solve_lower_cppad(%A%, %x%, %b%)%$$

$head Prototype$$
$srcfile%utility.hpp%
0%// BEGIN_SOLVE_LOWER_CPPAD_PROTOTYPE%// END_SOLVE_LOWER_CPPAD_PROTOTYPE%
1%$$

$head A$$
Is the a square lower triangular matrix in row-major order.
We use $icode n$$ to denote the number for rows and columns in $icode A$$;
i.e., $icode n$$ is equal to both $icode%A%.nr()%$$ and $icode%A%.nc()%$$.

$head x$$
This vector has size $icode n$$ and satisfies the equation
$codei%
	%A% * %x% = %b%
%$$

$head b$$
This is a vector with size $icode n$$.

$head msg$$
If the return value $icode msg$$ is empty, no error has occurred.
Otherwise, it is one of the following
$codei%
	"A is not in row-major order"
	"A is not lower triangular"
	"A is not invertible"
%$$.

$head Example$$
$srcfile%utility.cpp%
	0%// BEGIN_TEST_SOLVE_LOWER_CPPAD%// END_TEST_SOLVE_LOWER_CPPAD%$$

$end
*/
// BEGIN_SOLVE_LOWER_CPPAD_PROTOTYPE
template <class Scalar> std::string
solve_lower_cppad(
	const CPPAD_SPARSE(Scalar)&  A ,
	VECTOR(Scalar)&              x ,
	const VECTOR(Scalar)&        b )
// END_SOLVE_LOWER_CPPAD_PROTOTYPE
{	size_t nr = A.nr();
	size_t nnz = A.nnz();
	assert( nr == A.nc() );
	assert( nr == size_t( b.size() ) );
	assert( nr == size_t( x.size() ) );
	assert( nnz > 0 );
	//
	const VECTOR(size_t)& row( A.row() );
	const VECTOR(size_t)& col( A.col() );
	const VECTOR(Scalar)& val( A.val() );
	//
	//
	// check row-major order
	std::string msg = "A is not in row-major order";
	size_t k = 0;
	{	size_t r = row[k];
		size_t c = col[k];
		for(k = 1; k < nnz; k++)
		{	if( row[k] == r )
			{	if( col[k] <= c )
					return msg;
			}
			else if( row[k] < r )
				return msg;
			r = row[k];
			c = col[k];
		}
	}
	//
	// check lower triangular
	msg = "A is not lower triangular";
	for(k = 0; k < nnz; k++)
	{	if( row[k] < col[k] )
			return msg;
	}
	//
	// check for invertible and solve for x
	k = 0;
	msg = "A is not invertible";
	for(size_t i = 0; i < nr; i++)
	{	// solve for x[i]
		Scalar sum = b[i];
		//
		// x[c] has been determines for c < i
		while( row[k] == i && col[k] < i )
		{	sum -= val[k] * x[ col[k] ];
			k++;
		}
		//
		// determine x[i]
		if( row[k] != i )
			return msg;
		if( val[k] == Scalar(0.0) )
			return msg;
		//
		x[i] = sum / val[k];
		//
		k++;
	}
	msg = "";
	return msg;
}
/*
-------------------------------------------------------------------------------
$begin jac_constraint$$
$spell
	Jacobian
	jac
	xy
	CppAD
	rcv
$$
$section Compute Jacobian of Implicit Function Constraints$$

$head Syntax$$
$codei%jac_constraint(%L_y%, %L_fun%, %xy%, %work%)%$$

$head Prototype$$
$srcfile%utility.hpp
	%0%// BEGIN_JAC_CONSTRAINT_PROTOTYPE%// END_JAC_CONSTRAINT_PROTOTYPE%
1%$$

$head Purpose$$
We are given a function $latex L : \B{R}^n \timex \B{R}^m \rightarrow \B{R}^m$$
and define the implicit function $latex Y : \B{R}^n rightarrow \B{R}^m$$ by
the constraint equation
$latex \[
	L[ x , Y(x) ] = 0
\]$$
This routine computes the sparse derivative of $latex L (x , y)$$
w.r.t $latex y$$.

$head L_y$$
This argument must initially be created with the empty constructor; i.e.,
$codei%
	CPPAD_SPARSE(Scalar) %L_y%
%$$
It's return value is a sparse matrix representation of the Jacobian
of the constraints with respect to the implicit variables; i.e.,
the partial of $latex L(x, y)$$ w.r.t $latex y$$.
The return value can be used to speed up subsequent calls when
the operation sequence in $icode L_fun$$ does not change.
Note that this is an $latex m \times m$$ matrix; i.e.,
the column indices are relative to $latex y$$ and not for $latex (x, y)$$.
The entries in $icode L_y$$ are in row-major order.

$head L_fun$$
The operation sequence for the $cref constraint$$ function
$latex L : \B{R}^n \times \B{R}^m \rightarrow \B{R}^m$$
is stored in $icode L_fun$$.

$head xy$$
This is the argument $latex (x, y)$$ at which we are evaluating the Jacobian.
It size is $codei%n%+%m%$$.

$head work$$
This is a work vector used to reduce the work.
If the operation sequence ion $icode L_fun$$ changes,
$icode L_y$$ and $icode work$$ should be reset to empty.
Otherwise, they should remain unchanged from their previous return value.

$head Example$$

$subhead Simple$$
$srcfile%utility.cpp%
	0%// BEGIN_TEST_JAC_CONSTRAINT%// END_TEST_JAC_CONSTRAINT%$$

$subhead Control Problem$$
$srcfile%utility.cpp%
	0%// BEGIN_TEST_CONTROL_JAC_CONSTRAINT%// END_TEST_CONTROL_JAC_CONSTRAINT%$$

$end
*/
// BEGIN_JAC_CONSTRAINT_PROTOTYPE
template <class Scalar>
void jac_constraint(
	CPPAD_SPARSE(Scalar)&                             L_y         ,
	CppAD::ADFun<Scalar>&                             L_fun       ,
	VECTOR(Scalar)&                                   xy          ,
	CppAD::sparse_jac_work&                           work        )
// END_JAC_CONSTRAINT_PROTOTYPE
{	typedef CppAD::sparse_rc<VECTOR(size_t)> sparsity;
	size_t m      = L_fun.Range();
	size_t n      = L_fun.Domain() - m;
	//
	sparsity pattern_jac;
	if( L_y.nc() == 0 )
	{	// This work is only done when initial L_y is empty (first call)
		assert( L_y.nr() == 0 );
		assert( L_y.nnz() == 0 );
		//
		size_t nr     = n + m;
		size_t nc     = n + m;
		size_t nnz_in = n + m ;
		sparsity pattern_identity(nr, nc, nnz_in);
		for(size_t k = 0; k < nnz_in; k++)
		{	size_t r = k;
			size_t c = k;
			pattern_identity.set(k, r, c);
		}
		//
		// set pattern_jac to sparsity pattern for L^{(1)} (x, y)
		bool transpose     = false;
		bool dependency    = false;
		bool internal_bool = false;
		L_fun.for_jac_sparsity(
			pattern_identity, transpose, dependency, internal_bool, pattern_jac
		);
		//
		// sparsity pattern for L_y (x, y)
		const VECTOR(size_t)& row_jac( pattern_jac.row() );
		const VECTOR(size_t)& col_jac( pattern_jac.col() );
		VECTOR(size_t) row_major = pattern_jac.row_major();
		size_t nnz_jac    = pattern_jac.nnz();
		size_t nnz_subset = 0; // count number entries in L_y (x, y);
		for(size_t k = 0; k < nnz_jac; k++)
		{	size_t c = col_jac[k];
			if( c >= n )
				++nnz_subset;
		}
		sparsity pattern_shifted(m,   m, nnz_subset);
		size_t k_subset = 0;
		for(size_t k = 0; k < nnz_jac; k++)
		{	size_t ell = row_major[k];
			size_t c   = col_jac[ell];
			if( c >= n )
			{	size_t r = row_jac[ell];
				pattern_shifted.set(k_subset, r, c-n);
				++k_subset;
			}
		}
		assert( k_subset == nnz_subset );
		//
		// Set pattern in L_y and L_y
		L_y = CPPAD_SPARSE(Scalar)(pattern_shifted);
	}
	// L_y_subset
	size_t nnz = L_y.nnz();
	sparsity pattern_subset(m, n+m, nnz);
	for(size_t k = 0; k < nnz; k++)
		pattern_subset.set(k, L_y.row()[k], L_y.col()[k] + n);
	CPPAD_SPARSE(Scalar) L_y_subset(pattern_subset);
	//
	std::string coloring = "cppad";
	L_fun.sparse_jac_rev(
		xy, L_y_subset, pattern_jac, coloring, work
	);
	for(size_t k = 0; k < nnz; k++)
		L_y.set(k, L_y_subset.val()[k]);
	return;
}

# endif
