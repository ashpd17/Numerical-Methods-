# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
using namespace std;
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
double pythag ( double a, double b );
int tql2 ( int n, double d[], double e[], double z[] )
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
        if ( 1000 <= j )
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

double r8_abs ( double x )

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
//**************************80

double r8_epsilon ( )

{
  static double value = 2.220446049250313E-016;

  return value;
}
//**************************80

double r8_max ( double x, double y )

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
//**************************80

double r8_min ( double x, double y )

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
//**************************80

double r8_sign ( double x )

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
void extrapolate(double E[],double E2[],int n,double res[])
{
	for(int i = 0; i < n; i++)
	{
		res[i] = E2[i] + E2[i] - E[i];
	}
}
double p_f(double x)
{
	double pi = 3.14159265359;
	double value = 2 + sin(2*pi*x);
	return value;
}
double q_f(double x)
{
	double value = -10;
	return value;
}
double w_f(double x)
{
	double value = 1 + sqrt(x);
	return value;
} 
void fdm(int n, double a, double b, double alp_1,double alp_2,double bet_1, double bet_2,double P[])
{
	double x[n+1],e[n], h = (b-a)/n;
	for(int i = 0; i <= n; i++)
	{
		x[i] = a + i*h;
	}
	for(int i = 0; i < n; i++)
	{
		e[i] = a + (i + 0.5)*h;
	}
	double Q[n],W[n],A[n],B[n];
	for(int i = 0; i < n; i++)
	{
		Q[i] = q_f(e[i]);
		W[i] = w_f(e[i]);
	}
	double ya = 2 - 4*alp_2/(2*alp_2 - alp_1*h), yb = 2 - 4*bet_2/(2*bet_2 + bet_1*h);
	for(int i = 0; i < n; i++)
	{
		if(i == 0)
		{
			P[i] = 1/(pow(h,2))*(ya*p_f(x[i]) + p_f(x[i+1]));
			A[i] = -1/(pow(h,2))*p_f(x[i+1]);
			B[i] = -1/(pow(h,2))*p_f(x[i+1]);
		}
		else if(i == n - 1)
		{
			P[i] = 1/(pow(h,2))*(yb*p_f(x[i]) + p_f(x[i+1]));
		}
		else
		{
			A[i] = -1/(pow(h,2))*p_f(x[i+1]);
			B[i] = -1/(pow(h,2))*p_f(x[i+1]);
			P[i] = 1/(pow(h,2))*(p_f(x[i]) + p_f(x[i+1]));
		}
		
	}
	for(int i = 0; i < n; i++)
	{
		P[i] = P[i] + Q[i];
	}
	double L[n+1];
	for(int i = 0; i < n; i++)
	{
		L[i] = 1/sqrt(W[i]);
	}
	for(int i = 0; i < n-1; i++)
	{
		A[i] = A[i]*L[i+1]*L[i];
		B[i] = B[i]*L[i]*L[i+1];
	}
	for(int i=0; i < n; i++)    
	{
		P[i] = P[i]*L[i]*L[i];
	}
    double z[(n)*(n)];
	for(int i = 0; i < (n)*(n); i++)
	{
		if(i%(n+1) == 0)
		{
			z[i] = 1;
		}
		else
		{
			z[i] = 0;
		}
	}
	double E[n];
	E[0] = 0;
	for(int i = 1; i < n; i++)
	{
		E[i] = A[i-1];
	}
	tql2(n,P,E,z);
}
int main()
{
	cout.precision(10);    
    cout.setf(ios::fixed);
	double P[1000],P2[1000];
	int n = 500;
	double a = 0,b=1;
	double alp_1 = 1, bet_1 = 5, alp_2 = 0, bet_2 = 1;
	double res[n];
	fdm(n,a,b,alp_1,alp_2,bet_1,bet_2,P);
	n = n/2;
	fdm(n,a,b,alp_1,alp_2,bet_1,bet_2,P2);
	extrapolate(P,P2,2*n,res);
	for(int i = 0; i < 8; i++)
	{
		cout<<"lambda "<< i+1 <<" = "<<P[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<P2[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<res[i]<<endl;
	}
    return 0;
}
