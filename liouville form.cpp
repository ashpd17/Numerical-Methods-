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