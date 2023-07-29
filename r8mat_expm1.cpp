# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <ctime>
# include <cstring>

using namespace std;

# include "r8mat_expm1.h"
# include "r8lib.h"

//****************************************************************************80

double *r8mat_expm1 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 December 2011
//
//  Author:
//
//    Cleve Moler, Charles Van Loan
//
//  Reference:
//
//    Cleve Moler, Charles VanLoan,
//    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
//    Twenty-Five Years Later,
//    SIAM Review,
//    Volume 45, Number 1, March 2003, pages 3-49.
//
//  Parameters:
//
//    Input, int N, the dimension of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double R8MAT_EXPM1[N*N], the estimate for exp(A).
//
{
  double *a2;
  double a_norm;
  double c;
  double *d;
  double *e;
  int ee;
  int k;
  const double one = 1.0;
  int p;
  const int q = 6;
  int s;
  double t;
  double *x;

  a2 = r8mat_copy_new ( n, n, a );

  a_norm = r8mat_norm_li ( n, n, a2 );

  ee = ( int ) ( r8_log_2 ( a_norm ) ) + 1;
  
  s = i4_max ( 0, ee + 1 );

  t = 1.0 / pow ( 2.0, s );

  r8mat_scale ( n, n, t, a2 );

  x = r8mat_copy_new ( n, n, a2 );

  c = 0.5;

  e = r8mat_identity_new ( n );

  r8mat_add ( n, n, one, e, c, a2, e );

  d = r8mat_identity_new ( n );

  r8mat_add ( n, n, one, d, -c, a2, d );

  p = 1;

  for ( k = 2; k <= q; k++ )
  {
    c = c * ( double ) ( q - k + 1 ) / ( double ) ( k * ( 2 * q - k + 1 ) );

    r8mat_mm ( n, n, n, a2, x, x );

    r8mat_add ( n, n, c, x, one, e, e );

    if ( p )
    {
      r8mat_add ( n, n, c, x, one, d, d );
    }
    else
    {
      r8mat_add ( n, n, -c, x, one, d, d );
    }

    p = !p;
  }
//
//  E -> inverse(D) * E
//
  r8mat_minvm ( n, n, d, e, e );
//
//  E -> E^(2*S)
//
  for ( k = 1; k <= s; k++ )
  {
    r8mat_mm ( n, n, n, e, e, e );
  }

  delete [] a2;
  delete [] d;
  delete [] x;

  return e;
}
