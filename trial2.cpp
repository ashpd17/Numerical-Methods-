#include <iostream>
#include<fstream>
#include <cmath>
using namespace std;
double p_f(double x)
{
	double w = 1;
	return w;
}
double q_f(double x)
{
	double w = 0;
	return w;
}
double w_f(double x)
{
	double w = pow(1+x,-2);
	return w;
} 
int main()
{
	cout.precision(10);    
    cout.setf(ios::fixed);
	int n = 50;
	double a = 0,b=1;
	double x[n+1],e[n+1], h = (b-a)/n;
	//discretisation
	for(int i = 0; i < n; i++)
	{
		x[i] = a + i*h;
	}
	for(int i = 0; i < n; i++)
	{
		e[i] = a + (i + 1 - 0.5)*h;
	}
	double Q[n+1],W[n+1], P[n+1],A[n+1],B[n+1];
	for(int i = 0; i < n; i++)
	{
		Q[i] = q_f(e[i]);
		W[i] = w_f(e[i]);
	}
	double alp_1 = 1, bet_1 = 1, alp_2 = 0, bet_2 = 0;
	double ya = 2 - 4*alp_2/(2*alp_2 - alp_1*h), yb = 2 - 4*bet_2/(2*bet_2 + bet_1*h);
	for(int i = 0; i< n; i++)
	{
		if(i == 0)
		{
			P[i] = 1/(pow(h,2))*(ya*p_f(x[i]) + p_f(x[i+1]));
			A[i] = -1/(pow(h,2))*p_f(x[i]);
			B[i] = -1/(pow(h,2))*p_f(x[i]);
		}
		else if(i == n)
		{
			P[i] = 1/(pow(h,2))*(yb*p_f(x[i]) + p_f(x[i-1]));
		}
		else
		{
			A[i] = -1/(pow(h,2))*p_f(x[i]);
			B[i] = -1/(pow(h,2))*p_f(x[i]);
			P[i] = 1/(pow(h,2))*(p_f(x[i]) + p_f(x[i-1]));
		}
	}
	for(int i = 0; i < n; i++)
	{
		P[i] = P[i] + Q[i];
	}
	//cholsky decomposition
	double L[n+1];
	for(int i = 0; i < n; i++)
	{
		L[i] = 1/sqrt(W[i]);
	}
	for(int i = 0; i < n-1; i++)
	{
		A[i] = A[i]*L[i+1];
		B[i] = B[i]*L[i];
	}
	for(int i=0; i < n; i++)    
	{
		P[i] = P[i]*L[i];
	}
	for(int i=0; i < n; i++)    
	{
		P[i] = P[i]*L[i];
	}
	for(int i = 0; i < n-1; i++)
	{
		A[i] = A[i]*L[i];
		B[i] = B[i]*L[i+1];
	}
}
