//put all the working code here
#include<iostream>
#include "eigen.cpp"
using namespace std;
int main(){
int n1=3;
double d1[n1]={3,3,3};
double e1[n1]={1,1,1};
double z1[9]={1,0,0,0,1,0,0,0,1};
    cout<<tql2(n1,d1,e1,z1);
    cout<<d1[0]<<" "<<d1[1]<<" "<<d1[2]<<endl;
    for(int i=0;i<9;i++){
      cout<<"  " <<z1[i]<<" ";

    }
    system("pause");
    return 0;
}