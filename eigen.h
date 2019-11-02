#ifndef _EIGEN_H
#define _EIGEN_H
double dist ( double a, double b );//Returns hypotenuse of sides a and b
double eps();
double signval ( double x );//returns sign of input
int tqlrat ( int n, double d[], double e2[] );//calculates eigen value by rational tql2 method
int tql2 ( int n, double d[], double e[], double z[] );//calculates eigen value
#endif //_EIGEN_H