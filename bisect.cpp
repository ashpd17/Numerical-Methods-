# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
double epsilon = 2.220446049250313E-016
double bisect(double c[], double b[], double beta,int n, double  m1,double m2, double eps, double epsilon)
{
	double h,xmin,xmax;
	int i;
	beta[0] = b[0] = 0;
	//xmin , xmax calculation
	xmin = c[n] - abs(a[n]);
	xmax =  c[n] + abs(a[n]);
	for(int i = n-1; i > 0; i--)
	{
		double h = abs(b[i]) + abs(b[i+1]);
		if(c[i] + h > xmax)
		{
			xmax = c[i] + h;
		}
		if(c[i] - h < xmin)
		{
			xmin = c[i] - h;
		}
	}
	if(xmin + xmax > 0)
	{
		eps2 = epsilon*xmax;
	}
	else
	{
		eps2 = -xmin;
	}
	if(eps1 <= 0)
	{
		eps1 = eps2;
	}
	eps2 = 0.5*eps1 + 7*eps2;
	//inner block
	int a,k;
	double q,x1,xu,x0,wu[m1,m2];
	x0 = xmax;
	for(int i = m1; i < m2; i++)
	{
		x[i] = xmax;
		wu[i] = xmin;
		z = 0;
		//kth eigenvalue
		for(k = m2, k < m1; i--)
		{
			xu = xmin;
			for(int i = k; i > m1; i--)
			{
				if(xu < wu[i])
				{
					xu = xu[i];
					if(x0 > x[k])
					{
						x0 = x[k];
					}
					x1 = (xu + x0)/2;
					while(x0 - xu > 2*epsilon*(abs(xu) + abs(x0)+eps1)
					{
						z = z + 1;
						//sturm sequence
						a = 0;
						q
					double q = 1;
					for(int i = 0; i < n; i++)
					{
						if(q!= 0)
						{
							q = c[i] - x1 - beta[i]/q;
						}
						else
						{
							q = c[i] - x1 - abs(b[i]/epsilon);
						}
					}
					if(q < 0)
					{
						a = a + 1;
					}
					if(a < k)
					{
						if(a < m1)
						{
							xu = wu[m1] = x1;
						}
						else
						{
							xu = wu[a+1] = x1;
							if(x[a] > x1)
							{
								x[a] = x1;
							}
						}
					else
					{
						x0 = x1;
					}
					}
				}
				x[k] = (x0 + xu)/2;
			}
		}
	}
}