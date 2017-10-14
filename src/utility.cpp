# include "utility.hpp"

// ---------------------------------------------------------------------------
// BEGIN_TEST_SPARSE_CPPAD2EIGEN
bool test_sparse_cppad2eigen(void)
{	bool ok = true;
	//
	// create sparse_cppad
	size_t nr  = 6;
	size_t nc  = nr;
	size_t nnz = nr  + nr - 1;
	CppAD::sparse_rc<VECTOR(size_t)> pattern_in(nr, nc, nnz);
	for(size_t k = 0; k < nr; k++)
		pattern_in.set(k, k, k);
	for(size_t k = 0; k < nr - 1; k++)
		pattern_in.set(nr + k, k+1, k);
	CPPAD_SPARSE(double) sparse_in(pattern_in);
	for(size_t k = 0; k < nnz; k++)
		sparse_in.set(k, double(k+1));
	//
	SparseMatrix<double> sparse_out;
	sparse_cppad2eigen(sparse_in, sparse_out);
	//
	// check the result
	const VECTOR(size_t) row_in( sparse_in.row() );
	const VECTOR(size_t) col_in( sparse_in.col() );
	const VECTOR(double) val_in( sparse_in.val() );
	//
	typedef typename SparseMatrix<double>::InnerIterator iterator;
	VECTOR(size_t) col_major = sparse_in.col_major();
	size_t k = 0;
	for(int c = 0; c < sparse_out.outerSize(); ++c)
	{	for(iterator itr(sparse_out, c); itr; ++itr)
		{	size_t ell = col_major[k];
			ok &= size_t(itr.row()) == row_in[ell];
			ok &= size_t(itr.col()) == col_in[ell];
			ok &= itr.value() == val_in[ell];
			++k;
		}
	}
	//
	return ok;
}
// END_TEST_SPARSE_CPPAD2EIGEN
// ----------------------------------------------------------------------------
// BEGIN_TEST_SOLVE_LOWER_CPPAD
bool test_solve_lower_cppad(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	//     [ 1 0 0 ]
	// A = [ 2 3 0 ]
	//     [ 4 5 6 ]
	size_t nr  = 3;
	size_t nc  = 3;
	size_t nnz = 6;
	CppAD::sparse_rc< VECTOR(size_t) > pattern(nr, nc, nnz);
	size_t k = 0;
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = 0; j <= i; j++)
			pattern.set(k++, i, j);
	}
	assert( k == nnz );
	CPPAD_SPARSE(double) A(pattern);
	for(k = 0; k < nnz; k++)
		A.set(k, double(k+1));
	//
	// b = [1, 8, 32]^T
	VECTOR(double) b(nr);
	b[0] = 1.0;
	b[1] = 8.0;
	b[2] = 32.0;
	//
	// solve A * x = b
	VECTOR(double) x(nr);
	std::string msg = solve_lower_cppad<double>(A, x, b);
	//
	// check result
	ok &= msg == "";
	ok &= CppAD::NearEqual(x[0], 1.0, eps99, eps99);
	ok &= CppAD::NearEqual(x[1], 2.0, eps99, eps99);
	ok &= CppAD::NearEqual(x[2], 3.0, eps99, eps99);
	//
	return ok;
}
// END_TEST_SOLVE_LOWER_CPPAD
// ----------------------------------------------------------------------------
// BEGIN_TEST_JAC_CONSTRAINT
bool test_jac_constraint(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();

	// record z = L(x, y)
	VECTOR( CppAD::AD<double> ) axy(2), az(1);
	axy[0] = 0.5;
	axy[1] = 0.5;
	CppAD::Independent( axy );
	az[0] = axy[0] * axy[0] + axy[1] * axy[1] - 1.0;
	CppAD::ADFun<double> L_fun(axy, az);

	// choose a point that satisfies the constraints
	VECTOR(double) xy(2);
	xy[0] = 0.5;
	xy[1] = std::sqrt( 1.0 - xy[0] * xy[0] );

	// Evaluate L_y (x, y)
	CPPAD_SPARSE(double) L_y;
	CppAD::sparse_jac_work work;
	jac_constraint(L_y, L_fun, xy, work);

	// check the result is L_y(x, y) = 2.0 * y
	ok &= L_y.nr()  == L_fun.Range();
	ok &= L_y.nc()  == L_fun.Range();
	ok &= L_y.nnz() == 1;
	ok &= L_y.row()[0] == 0;
	ok &= L_y.col()[0] == 0;
	ok &= CppAD::NearEqual( L_y.val()[0], 2.0 * xy[1], eps99, eps99 );
	//
	return ok;
}
// END_TEST_JAC_CONSTRAINT
// ----------------------------------------------------------------------------
// BEGIN_TEST_CONTROL_JAC_CONSTRAINT
# include "control.hpp"
bool test_control_jac_constraint(void)
{	bool ok = true;
	double eps99 = 99.0 * std::numeric_limits<double>::epsilon();
	//
	// record L_fun
	CppAD::ADFun<double>       L_fun;
	size_t                     J = 6;
	double                     delta_t = 2.0;
	VECTOR(double) p(4);
	for(size_t i = 0; i < 4; i++)
		p[i] = double(5 - i);
	control::rec_constraint(L_fun, J, delta_t, p);
	size_t m      = L_fun.Range();
	size_t n      = L_fun.Domain() - m;
	//
	// compute L_y (x, y)
	CPPAD_SPARSE(double) L_y;
	VECTOR(double) xy(n + m);
	for(size_t i = 0; i < n+m; i++)
		xy[i] = 0.5;
	CppAD::sparse_jac_work work;
	//
	// all vj are 0.5
	double vj  = 0.5;
	double r   = sqrt( (vj + 1.0)*(vj + 1.0) + vj * vj );
	double r3  = r * r * r;
	double r5  = r * r * r * r * r;
	double a   = 1.0 / r3 - 1.0;
	double ap0 = -3.0 * (vj + 1.0) / r5;
	double ap1 = -3.0 * vj / r5;
	//
	// check repeated calls to jac_constraint
	for(size_t count = 0; count < 2; count ++)
	{	jac_constraint(L_y, L_fun, xy, work);
		//
		// check Jacobian
		size_t nnz = L_y.nnz();
		const VECTOR(size_t)& row( L_y.row() );
		const VECTOR(size_t)& col( L_y.col() );
		const VECTOR(double)& val( L_y.val() );
		for(size_t k = 0; k < nnz; k++)
		{	//
			size_t i = row[k];
			size_t j = col[k];
			double v = val[k];
			//
			// check lower triangular
			ok &= j <= i;
			//
			bool diagonal_block = (i / 4) == (j / 4);
			if( diagonal_block )
			{	ok &= i == j;
				ok &= v == 1.0;
			}
			else
			{	double check;
				// must be the lower diagonal block
				ok &= (i / 4) == (j / 4) + 1;
				//
				// derivative of this v w.r.t previous
				size_t i_4 = i % 4;
				size_t j_4 = j % 4;
				if( i_4 == 0 )
				{	// derivative of this v0 w.r.t previous v0 or v2
					ok &= j_4 == 0 || j_4 == 2;
					if( j_4 == i_4 )
						ok &= v == -1.0;
					else
						ok &= v == - 1.0 * delta_t;
				}
				if( i_4 == 1 )
				{	// derivative of this v1 w.r.t previous v1 or v3
					ok &= j_4 == 1 || j_4 == 3;
					if( j_4 == i_4 )
						ok &= v == -1.0;
					else
						ok &= v == - 1.0 * delta_t;
				}
				if( i_4 == 2 )
				{	// derivative of this v2
					switch( j_4 )
					{
						case 0:
						// w.r.t previous v0
						check = (a + (1.0 + vj) * ap0) * delta_t;
						ok &= CppAD::NearEqual(v, check, eps99, eps99);
						break;

						case 1:
						// w.r.t previous v1
						check = (1.0 + vj) * ap1 * delta_t;
						ok &= CppAD::NearEqual(v, check, eps99, eps99);
						break;

						case 2:
						// w.r.t previous v2
						ok &= v == -1.0;
						break;

						case 3:
						// w.r.t previous v3
						ok &= v == - 2.0 * delta_t;
						break;

						default:
						ok = false;
					}
				}
				if( i_4 == 3 )
				{	// derivative of this v3
					switch( j_4 )
					{
						case 0:
						// w.r.t previous v0
						check = vj * ap0 * delta_t;
						ok &= CppAD::NearEqual(v, check, eps99, eps99);
						break;

						case 1:
						// w.r.t previous v1
						check = (a + vj * ap1) * delta_t;
						ok &= CppAD::NearEqual(v, check, eps99, eps99);
						break;

						case 2:
						// w.r.t previous v2
						ok &= v == 2.0 * delta_t;
						break;

						case 3:
						// w.r.t previous v3
						ok &= v == -1.0;
						break;

						default:
						ok = false;
					}
				}
			}
		}
	}
	return ok;
}
// END_TEST_CONTROL_JAC_CONSTRAINT

