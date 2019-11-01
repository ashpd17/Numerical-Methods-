//including header files
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include"eigen.cpp"//contains tql2 method
#include"companion.cpp"//contains other functions
using namespace std;
int main(){
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
	extrapolate(P,P2,2*n,res);//for richardson extrapolation
	for(int i = 0; i < 8; i++){
		cout<<"lambda "<< i+1 <<" = "<<P[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<P2[i]<<"\t";
		cout<<"lambda "<< i+1 <<" = "<<res[i]<<endl;
	}
    system("pause");
    return 0;
}