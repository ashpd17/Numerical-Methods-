#include"eigen.h"
using namespace std;
int tqlrat ( int n, double d[], double e2[] )
{
  double b = 0,c = 0,f,g,h;
  int i, errcd, ii, j, l, l1, m, mml;
  double p, r, s, t;
  errcd = 0;
  if (n == 1){
      return errcd;
    }
  for ( i = 1; i < n; i++ ){
    e2[i-1] = e2[i];
    }

  f = 0.0;
  t = 0.0;
  e2[n-1] = 0.0;
  for ( l = 0; l < n; l++ ){
     j = 0;
     h = abs ( d[l] ) + sqrt ( e2[l] );
     if ( t <= h ){
       t = h;
       b = abs ( t ) * eps ( );
       c = b * b;
    }
//  Look for small squared sub-diagonal element
for ( m = l; m < n; m++ ){  
      if ( e2[m] <= c )
      {
        break;
      }
    }

    if ( m != l ){
        //edit01
      while(true)      {
        if ( 50 <= j ){
          errcd = l + 1;
          return errcd;
        }

        j++;
        //shifting value position
        l1 = l + 1;
        s = sqrt ( e2[l] );
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * s );
        r = dist ( p, 1.0 );
        d[l] = s / ( p + abs ( r ) * signval ( p ) );
        h = g - d[l];
        for ( i = l1; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
//  Rational ql subroutine.

        g = d[m];
        if ( g == 0.0 )
        {
          g = b;
        }

        h = g;
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          i = m - ii;
          p = g * h;
          r = p + e2[i];
          e2[i+1] = s * r;
          s = e2[i] / r;
          d[i+1] = h + s * ( h + d[i] );
          g = d[i] - e2[i] / g;
          if ( g == 0.0 )
          {
            g = b;
          }
          h = g * p / r;
        }
        e2[l] = s * g;
        d[l] = h;
//cheking if value gets zero
        if ( h == 0.0 )
        {
          break;
        }

        if ( abs ( e2[l] ) <= abs ( c / h ) )
        {
          break;
        }

        e2[l] = h * e2[l];

        if ( e2[l] == 0.0 )
        {
          break;
        }
      }
    }

    p = d[l] + f;
//Ordering the vector 
    for ( i = l; 0 <= i; i-- )
    {
      if ( i == 0 )
      {
        d[i] = p;
        break;
      }
      else if ( d[i-1] <= p )
      {
        d[i] = p;
        break;
      }
      d[i] = d[i-1];
    }
  }

  return errcd;
}

int tql2 ( int n, double d[], double e[], double z[] ){
  double c, c2, c3 = 0, dl1, el1, f, g, h;
  int i, errcd, ii, j, k, l, l1, l2, m, mml;
  double p, r, s, s2 = 0, t, tst1, tst2;
  errcd = 0;
  if ( n == 1 ){
    return errcd;
  }

  for ( i = 1; i < n; i++ ){
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  
  for ( l = 0; l < n; l++ ){
    j = 0;
    h = abs ( d[l] ) + abs ( e[l] );
    tst1 = max ( tst1, h );
//Finding smallest sub diagonal element
    for ( m = l; m < n; m++ ){
      tst2 = tst1 + abs ( e[m] );
      if ( tst2 == tst1 ){
        break;
      }
    }

    if ( m != l ){
      for ( ; ; ){
        if ( 1000 <= j ){
          errcd = l + 1;
          return errcd;
        }

        j = j + 1;
//shifting by 1 position
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = dist ( p, 1.0 );
        d[l] = e[l] / ( p + signval ( p ) * abs ( r ) );
        d[l1] = e[l] * ( p + signval ( p ) * abs ( r ) );
        dl1 = d[l1];
        h = g - d[l];
        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
//  QL transformation.
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ ){
          c3 = c2;
          c2 = c;
          s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = dist ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
//Initializing vectors
          for ( k = 0; k < n; k++ ){
            h = z[k+(i+1)*n];
            z[k+(i+1)*n] = s * z[k+i*n] + c * h;
            z[k+i*n] = c * z[k+i*n] - s * h;
          }
        }
        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + abs ( e[l] );

        if ( tst2 <= tst1 ){
          break;
        }
      }
    }
    d[l] = d[l] + f;
  }
//Ordering eigen values and eigen vectors
  for ( ii = 1; ii < n; ii++ ){
    i = ii - 1;
    k = i;
    p = d[i];
    for ( j = ii; j < n; j++ ){
      if ( d[j] < p ){
        k = j;
        p = d[j];
      }
    }

    if ( k != i ){
      d[k] = d[i];
      d[i] = p;
      for ( j = 0; j < n; j++ ){
        t        = z[j+i*n];
        z[j+i*n] = z[j+k*n];
        z[j+k*n] = t;
      }
    }
  }
  return errcd;
}
double dist ( double a, double b ){//Returns hypotenuse of right triangle with side a and b
  double p, r, s, t, u;
  p = max ( abs ( a ), abs ( b ) );//choosing larger value
  if ( p != 0.0 ){
    r = min ( abs ( a ), abs ( b ) ) / p;
    r = r * r;
    while ( 1 ){
      t = 4.0 + r;
      if ( t == 4.0 ){
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
double eps ( ){
  static double value = 2.220446049250313E-016;
  return value;
}
double signval ( double x ){
  double value;
 ( x < 0.0 )?value = -1.0:value = 1.0;
  return value;
}