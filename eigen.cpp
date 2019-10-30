# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include "eigen.h"
using namespace std;
int tql2 ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues and eigenvectors of a symmetric
//    tridiagonal matrix by the QL method.  The eigenvectors of a full
//    symmetric matrix can also be found if TRED2 has been used to reduce this
//    full matrix to tridiagonal form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Bowdler, Martin, Reinsch, Wilkinson,
//    TQL2,
//    Numerische Mathematik,
//    Volume 11, pages 293-306, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, the diagonal elements of
//    the matrix.  On output, the eigenvalues in ascending order.  If an error
//    exit is made, the eigenvalues are correct but unordered for indices
//    1,2,...,IERR-1.
//
//    Input/output, double E[N].  On input, E(2:N) contains the
//    subdiagonal elements of the input matrix, and E(1) is arbitrary.
//    On output, E has been destroyed.
//
//    Input, double Z[N*N].  On input, the transformation matrix
//    produced in the reduction by TRED2, if performed.  If the eigenvectors of
//    the tridiagonal matrix are desired, Z must contain the identity matrix.
//    On output, Z contains the orthonormal eigenvectors of the symmetric
//    tridiagonal (or full) matrix.  If an error exit is made, Z contains
//    the eigenvectors associated with the stored eigenvalues.
//
//    Output, int TQL2, error flag.
//    0, normal return,
//    J, if the J-th eigenvalue has not been determined after
//    30 iterations.
//
{
  double c;
  double c2;
  double c3 = 0;
  double dl1;
  double el1;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int k;
  int l;
  int l1;
  int l2;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double s2 = 0;
  double t;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    j = 0;
    h = r8_abs ( d[l] ) + r8_abs ( e[l] );
    tst1 = r8_max ( tst1, h );
//
//  Look for a small sub-diagonal element.
//
    for ( m = l; m < n; m++ )
    {
      tst2 = tst1 + r8_abs ( e[m] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
//
//  Form shift.
//
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = pythag ( p, 1.0 );
        d[l] = e[l] / ( p + r8_sign ( p ) * r8_abs ( r ) );
        d[l1] = e[l] * ( p + r8_sign ( p ) * r8_abs ( r ) );
        dl1 = d[l1];
        h = g - d[l];
        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
//
//  QL transformation.
//
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = pythag ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
//
//  Form vector.
//
          for ( k = 0; k < n; k++ )
          {
            h = z[k+(i+1)*n];
            z[k+(i+1)*n] = s * z[k+i*n] + c * h;
            z[k+i*n] = c * z[k+i*n] - s * h;
          }
        }
        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + r8_abs ( e[l] );

        if ( tst2 <= tst1 )
        {
          break;
        }
      }
    }
    d[l] = d[l] + f;
  }
//
//  Order eigenvalues and eigenvectors.
//
  for ( ii = 1; ii < n; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i];
    for ( j = ii; j < n; j++ )
    {
      if ( d[j] < p )
      {
        k = j;
        p = d[j];
      }
    }

    if ( k != i )
    {
      d[k] = d[i];
      d[i] = p;
      for ( j = 0; j < n; j++ )
      {
        t        = z[j+i*n];
        z[j+i*n] = z[j+k*n];
        z[j+k*n] = t;
      }
    }
  }
  return ierr;
}
double pythag ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PYTHAG computes SQRT ( A * A + B * B ) carefully.
//
//  Discussion:
//
//    The formula
//
//      PYTHAG = sqrt ( A * A + B * B )
//
//    is reasonably accurate, but can fail if, for example, A^2 is larger
//    than the machine overflow.  The formula can lose most of its accuracy
//    if the sum of the squares is very large or very small.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Modified:
//
//    08 November 2012
//
//  Parameters:
//
//    Input, double A, B, the two legs of a right triangle.
//
//    Output, double PYTHAG, the length of the hypotenuse.
//
{
  double p;
  double r;
  double s;
  double t;
  double u;

  p = r8_max ( r8_abs ( a ), r8_abs ( b ) );

  if ( p != 0.0 )
  {
    r = r8_min ( r8_abs ( a ), r8_abs ( b ) ) / p;
    r = r * r;

    while ( 1 )
    {
      t = 4.0 + r;

      if ( t == 4.0 )
      {
        break;
      }

      s = r / t;
      u = 1.0 + 2.0 * s;
      p = u * p;
      r = ( s / u ) * ( s / u ) * r;
    }
  }
  return p;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  static double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
