//including header files
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <fstream>
#include"eigen.cpp"//contains tql2 method
#include"companion.cpp"//contains other functions
using namespace std;
double p_f(double x)
{
		return 1.0;
}
double q_f(double x)
{
	double value =(0.25);
	return 0;
}
double w_f(double x)
{   double pi = 3.14159265359;
	double value = (64*pow(pi,2)/(9*pow(x,6)));
	return 1.0;
}
ofstream ofs("lambda.csv");
int main(){
	cout.precision(10);    
    cout.setf(ios::fixed);
	double P[1000],P2[1000];
	int n = 200;
	double a = 0,b=1.0;
	double alp_1 = 1, alp_2 = 0, bet_1 = 1, bet_2 = 0;
	double res[n];
  	double dk[n*n];
    double dk1[n*n];
	for(int i = 0; i < (n)*(n); i++)
	{
		if(i%(n+1) == 0)
		{
			dk[i] = 1;
            dk1[i] = 1;
		}
		else
		{
			dk[i] = 0;
            dk1[i] = 0;
		}
	}
	fdm(n,a,b,alp_1,alp_2,bet_1,bet_2,P,dk,p_f,q_f,w_f);
	n = n/2;
	fdm(n,a,b,alp_1,alp_2,bet_1,bet_2,P2,dk1,p_f,q_f,w_f);
	extrapolate(P,P2,2*n,res);//for richardson extrapolation
	for(int i = 0; i < 8; i++){
		cout<<"lambda "<< i+1 <<" = "<<P[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<P2[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<res[i]<<"\t";
		cout<<"lambda"<<i+1<<" = "<<pow((i+1)*3.14,2)<<endl;
	}
	n=n*2;
    for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
			//cout<<j/(100.0)<<","<<dk[i*n+j]<<endl;
			ofs<<j/(200.0)<<","<<dk[i*n+j]<<endl;
		}
        
    }
	ofs.close();
    system("pause");
    return 0;
}