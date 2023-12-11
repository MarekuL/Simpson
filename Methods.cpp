#include"Simpson.h"
double  F1(double x)
{
	return(sqrt(1 + 2*pow(x,3)));
}
double F2(double x, double y)
{
	return(4 - pow(x, 2) - pow(y, 2));
}
double Trapez(double a, double b)
{
	double e1 = 1E-4;
	double h = (b - a) / 2;
	double f0 = F1(a);
	double fn = F1(b);
	double sum1 = 0;
	double sum2 = 0;
	int k = 1;
	for (k;;)
	{
		sum2 = sum1;
		int n = (b - a) / h;
		for (int i = 1; i < n - 1; i++)
		{
			sum1 += 2 * F1(a + h * i);
		}
		sum1 += f0;
		sum1 += fn;
		sum1 *= (h / 2);
		h /= 2;
		k++;
		if (abs(sum1 - sum2) <= 3 * e1)
			break;
	}
	cout << "Trapez: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << e1;

	return sum1;
}
double Simpson(double a, double b)
{

	double e2 = 1E-5;
	double h = (b - a) / 2;
	double f0 = F1(a);
	double fn = F1(b);
	double sum1 = 0;
	double sum2 = 0;
	int k = 1;
	for (k;;)
	{
		sum2 = sum1;

		double n = (b - a) / h;

		sum1 += F1(a);
		sum1 += F1(b);
		for (int i = 2; i < n - 1; i += 2)
		{
			sum1 += 2 * F1(a + h * i);
		}

		for (int i = 1; i < n; i += 2)
		{
			sum1 += 4 * F1(a + h * (i - 1));
		}
		sum1 /= 3;
		sum1 *= h;


		h /= 2;
		k++;
		if (abs(sum1 - sum2) <= 15 * e2)
			break;
	}
	cout << "Simpson: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << e2;

	return sum1;
}
double CubeSimpson(double a, double b, double c, double d)
{
	double hx = (b - a) / 2;
	double hy = (d - c) / 2;
	double sum1 = 0;
	double sum2 = 0;
	double e2 = 1E-5;
	int k = 0;



	for (k;;)
	{

		sum2 = sum1;
		sum1 = 0;

		double n = (b - a) / hx;
		double m = (d - c) / hy;

		for (int i = 0; i < n - 1; i += 2)
		{
			for (int j = 0; j < m - 1; j += 2)
			{
				sum1 += F2(a + hx * i, c + hy * j);
				sum1 += 4 * F2(a + hx * (i + 1), c + hy * j);
				sum1 += F2(a + hx * (i + 2), c + hy * j);
				sum1 += 4 * F2(a + hx * i, c + hy * (j + 1));
				sum1 += 16 * F2(a + hx * (i + 1), c + hy * (j + 1));
				sum1 += 4 * F2(a + hx * (i + 2), c + hy * (j + 1));
				sum1 += F2(a + hx * i, c + hy * (j + 2));
				sum1 += 4 * F2(a + hx * (i + 1), c + hy * (j + 2));
				sum1 += F2(a + hx * (i + 2), c + hy * (j + 2));
			}
		}

		sum1 *= hx * hy;
		sum1 /= 9;

		hx /= 2;
		hy /= 2;

		k++;



		if (abs(sum1 - sum2) <= 15 * e2)
			break;

	}

	cout << "SimpsonCube: " << sum1 << " K: " << k << " Del: " << sum1 - sum2 << " e: " << e2;

	return sum1;
}