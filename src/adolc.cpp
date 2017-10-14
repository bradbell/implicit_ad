
# include <adolc/adolc.h>
int main(void)
{
	int n = 1;    // number of independent and dependent variables
	int d = 2;    // highest order derivative
	int p = 1;    // number of directions
	int size = 3; // [(p + d) choose d] = (p+d-0) *...* (p+1) / (1 *...* d)
	//
	double* S[1];      // seed matrix S[n][p]
	S[0] = new double[p];
	//
	double* tensor[1]; // partials of z w.r.t x tensor[n][size]
	tensor[0] = new double[size];

	// Set Seed matrix to identity
	S[0][0] = 1.0;

	// independent and dependent variables
	adouble az[1], ay[1]; // az[n], ay[n]
	double  z[1],  y[1];  //  z[n],  y[n]

	// record operations for y = F(z) = sin(z)
	short int tag = 1;     // tape identifier
	trace_on(tag);         // start recording
	z[0]  = 1.0;           // value during recording
	az[0] <<= z[0];        // independent variable
	ay[0] = sin(az[0]);    // function evaluation
	ay[0] >>= y[0];        // dependent variable
	trace_off();           // turn off recording

	// evaluate inverse of F at x0
	double x0 = 0.5;            // argument for inverse
	z[0]      = std::asin(x0);  // F^-1 (x0) = z0
	inverse_tensor_eval(tag, n, d, p, z, tensor, S);

	// print tensor
	for(int j = 0; j < size; j++)
	{	std::cout << "tensor[0][" << j << "] = " << tensor[0][j] << "\n";
	}

	// derivative of F^-1 (x0)  =  1.0 / sqrt( 1 - x0 * x0 )
	double zp = 1.0 / std::sqrt(1.0 - x0 * x0 );

	// second derivative of F^-1 (x0)  =  x0 / sqrt( 1 - x0 * x0 )^3
	double zpp = x0 * zp * zp * zp;

	// print derivatives for z
	std::cout << "z = " << z[0] << "\n";
	std::cout << "zp = " << zp << "\n";
	std::cout << "zpp = " << zpp << "\n";
	//
	delete [] S[0];
	delete [] tensor[0];
	//
	std::cout << "adolc: Done\n";
	return 0;
}
