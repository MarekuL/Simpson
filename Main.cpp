#include"Simpson.h"
int main()
{


	double a = 1.2;
	double b = 2.471, c, d;
	double h = (b - a) / 2;
	double f0 = F1(a);
	double fn = F1(b);
	double sum1 = 0;
	double sum2 = 0;
	cout << endl;
	Trapez(a, b);
	cout << endl << endl;
	Simpson(a, b);
	cout << endl << endl;
	a = -1.0; b = 1.0; c = -1.0; d = 1.0;
	CubeSimpson(a, b, c, d);
	cout << endl;



}