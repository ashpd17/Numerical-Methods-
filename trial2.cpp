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
	double w = 1;
	return w;
} 
int main()
{
	cout.precision(10);    
    cout.setf(ios::fixed);
	int n = 5;
	double a = 0,b=1;
	double x[n+1],e[n+1], h = (b-a)/n;
	for(int i = 0; i < n; i++)
	{
		x[i] = a + i*h;
	}
	for(int i = 0; i < n; i++)
	{
		e[i] = a + (i + 1 - 0.5)*h;
		cout<<e[i]<<endl;
	}
	double Q[n+1][n+1],W[n+1][n+1], P[n+1][n+1],A[n+1][n+1];
	for(int i = 0; i< n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(i==j)
			{
				Q[i][j] = q_f(e[i]);
				W[i][j] = w_f(e[i]);
			}
			else
			{
				Q[i][j] = 0;
				W[i][j] = 0;
			}
		}
	}
	double alp_1 = 1, bet_1 = 1, alp_2 = 0, bet_2 = 0;
	double ya = 2 - 4*alp_2/(2*alp_2 - alp_1*h), yb = 2 - 4*bet_2/(2*bet_2 + bet_1*h);
	for(int i = 0; i< n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(i == j && i == 0)
			{
				P[i][j] = 1/(pow(h,2))*(ya*p_f(x[i]) + p_f(x[i+1]));
			}
			else if(i == j && i == n-1)
			{
				P[i][j] = 1/(pow(h,2))*(yb*p_f(x[i]) + p_f(x[i-1]));
			}
			else if(i == j && i != 0 && i != n)
			{
				P[i][j] = 1/(pow(h,2))*(p_f(x[i]) + p_f(x[i-1]));
			}
			else if(abs(i - j) == 1)
			{
				P[i][j] = -1/(pow(h,2))*p_f(x[i]);
			}
			else
			{
				P[i][j] = 0;
			}
		}
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			A[i][j] = P[i][j] + Q[i][j];
		}
	}
	for(int i = 0; i< n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout<<A[i][j]<<"\t";
		}
		cout<<endl;
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			if(i == j)
			{
				W[i][j] = 1/W[i][j];
			}
		}
	}
	double mul[n+1][n+1];
	for(int i=0; i<n; i++)    
	{    
		for(int j=0;j<n;j++)    
		{    
			mul[i][j]=0;    
			for(int k=0;k<n;k++)    
			{    
				mul[i][j]+=A[i][k]*W[k][j];    
			}    
		}	    
	}
	for(int i = 0; i< n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout<<mul[i][j]<<"\t";
		}
		cout<<endl;
	}
}