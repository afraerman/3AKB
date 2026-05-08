#pragma once
#include "alltypes.h"
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


#define IND 960
#define BIG std::pow(2.0, IND)
#define BIGI std::pow(2.0, -1*IND)
#define BIGS std::pow(2.0, IND / 2)
#define BIGIS std::pow(2.0, -1*IND / 2)
#define root3 1.7320508075688773


struct X{
	double x;
	int ix;
};

void xnorm(X& x)
{
	double w = std::fabs(x.x);
	if (w >= BIGS)
	{
		x.x = x.x * BIGI;
		x.ix++;
	}
	else if (w < BIGIS)
	{
		x.x = x.x * BIG;
		x.ix--;
	}
}

double x2f(X x)
{
	if (x.ix == 0) return x.x;
	if (x.ix < 0) return x.x * BIGI;
	return x.x * BIG;
}

X mul(double f, X x)
{
	X y;
	y.x = x.x * f;
	y.ix = x.ix;
	xnorm(y);
	return y;
}

X add(double f, X x, double g, X y)
{
	X z;
	int id = x.ix - y.ix;
	if (id == 0)
	{
		z.x = f * x.x + g * y.x;
		z.ix = x.ix;
	}
	else if (id == 1)
	{
		z.x = f * x.x + g * (y.x * BIGI);
		z.ix = x.ix;
	}
	else if (id == -1)
	{
		z.x = g * y.x + f * (x.x * BIGI);
		z.x = y.ix;
	}
	else if (id > 1)
	{
		z.x = f * x.x;
		z.ix = x.ix;
	}
	else
	{
		z.x = g * y.x;
		z.ix = y.ix;
	}
	xnorm(z);
	return z;
}

void all_pnm(double u, double t, int N, double cl, double sl, double** P, double** P1) // u = cos(phi), t = sin(phi)
{
	double d, a, b, f;
	int n, m;
	// store all pnm as f-numbers
	//double** P = new double*[N+1];
	//for (int i = 0; i < N+1; i++) P[i] = new double[N+1]();

	int M = N;

	P[0][0] = 1.0;
	P[1][0] = root3 * t;
	P[1][1] = root3 * u;

	P1[0][0] = 0.0;
	P1[1][0] = (t * P[1][0] - root3 * P[0][0]) / u;
	P1[1][1] = t * P[1][1] / u;

	// for P[n][0] no need to use X-numbers
	m = 0;
	for (n = m + 2; n < N+1; n++)
	{
		a = std::sqrt((double)((2.0 * n + 1.0) * (2.0 * n - 1)) / (double)((n + m) * (n - m)));
		b = -1.0 * std::sqrt((double)((2.0*n+1)*(n+m-1)*(n-m-1)) / (double)((2.0*n-3)*(n-m)*(n+m)));
		f = -1.0 * std::sqrt((double)((n-m)*(n+m)*(2*n+1)) / (double)(2*n-1));
		P[n][0] = a * t * P[n-1][0] + b * P[n-2][0];
		P1[n][0] = ((double)n * t * P[n][0] + f * P[n-1][0]) / u;
	}
	P[N+1][0] = 1.0; // cos(0 * lambda)
	P[N+2][0] = 0.0; // sin(0 * lambda)

	// Create all pmm as X-numbers
	X* pmm = new X[M+1];
	X x;
	x.x = root3 * u;
	x.ix = 0;

	X y{0.0, 0}, z{0.0, 0}, z1{0.0, 0};

	pmm[1] = x;
	// P[n][1] loop
	m = 1;
	y = x;
	n = m + 1;
	a = std::sqrt((double)(2.0 * n + 1));
	f = -1.0 * std::sqrt((double)((n-m)*(n+m)*(2*n+1)) / (double)(2*n-1));
	x = mul(a * t, y);
	xnorm(x);
	P[n][m] = x2f(x);
	P1[n][m] = ((double)n * t * P[n][m] + f * P[n-1][m]) / u;
	for (n = m + 2; n < N+1; n++)
	{
		a = t * std::sqrt((double)((2.0 * n + 1.0) * (2.0 * n - 1)) / (double)((n + m) * (n - m)));
		b = -1.0 * std::sqrt((double)((2.0*n+1)*(n+m-1)*(n-m-1)) / (double)((2.0*n-3)*(n-m)*(n+m)));
		f = -1.0 * std::sqrt((double)((n-m)*(n+m)*(2*n+1)) / (double)(2*n-1));
		z = add(a, x, b, y);
		P[n][m] = x2f(z);
		z1 = add((double)n * t / u, z, f / u, x);
		P1[n][m] = x2f(z1);
		y = x;
		x = z;
	}
	P[N+1][1] = cl;
	P[N+2][1] = sl;

	x = pmm[m];

	// P[m][m] as X-numbers
	for (m = 2; m < M+1; m++)
	{
		P[N+1][m] = P[N+1][m-1] * cl - P[N+2][m-1] * sl; // cos(m * lambda) = cos((m-1)*lambda) * cos(lambda) - sin((m-1)*lambda) * sin(lambda)
		P[N+2][m] = P[N+2][m-1] * cl + P[N+1][m-1] * sl; // sin(m * lambda) = sin((m-1)*lambda) * cos(lambda) + cos((m-1)*lambda) * sin(lambda)

		d = std::sqrt((double)(2.0 * m + 1) / (double)(2.0 * m));
		x = mul(d * u, x);
		xnorm(x);
		pmm[m] = x;

		P[m][m] = x2f(x);
		P1[m][m] = (double)m * t / u * P[m][m];
		if (m >= N) break;

		y = x;
		n = m + 1;
		a = std::sqrt((double)(2.0 * n + 1));
		b = 0.0;
		f = -1.0 * std::sqrt((double)(2.0 * m + 3.0));
		x = mul(a * t, y);
		xnorm(x);
		P[n][m] = x2f(x);
		P1[n][m] = ((double)n * t * P[n][m] + f * P[n-1][m]) / u;

		// P[n][m] as F-numbers
		for (n = m + 2; n < N+1; n++)
		{
			a = std::sqrt((double)((2.0 * n + 1.0) * (2.0 * n - 1)) / (double)((n + m) * (n - m)));
			b = -1.0 * std::sqrt((double)((2.0*n+1)*(n+m-1)*(n-m-1)) / (double)((2.0*n-3)*(n-m)*(n+m)));
			f = -1.0 * std::sqrt((double)((n-m)*(n+m)*(2*n+1)) / (double)(2*n-1));
			z = add(a * t, x, b, y);
			P[n][m] = x2f(z);
			z1 = add((double)n * t / u, z, f / u, x);
			P1[n][m] = x2f(z1);
			y = x;
			x = z;
		}

		x = pmm[m];
	}
	return;
}

void gravityFukushima(double r, double lat, double lon, int nmax, std::array<double, 3>& Result)
{
	double lat_rad = lat * M_PI / 180.0;
	double lon_rad = lon * M_PI / 180.0;

	double u = std::cos(lat_rad);
	double t = std::sin(lat_rad);

	double cl = std::cos(lon_rad);
	double sl = std::sin(lon_rad);

	double** P = new double* [nmax+3];
	for (int i = 0; i < nmax+3; i++) P[i] = new double[nmax+1];

	double** P1 = new double* [nmax+1];
	for (int i = 0; i < nmax+1; i++) P1[i] = new double[nmax+1];

	all_pnm(u, t, nmax, cl, sl, P, P1);
	
	/*
	std::cout << "\033[31m";
	for (int n = 0; n <= nmax; n++)
	{
		for (int m = 0; m <= n; m++)
		{
			std::cout << P[n][m] << "\t|";
		}
		std::cout << std::endl;
	}
	std::cout << "\033[0m";
	*/

	/*
	std::cout << "\033[31m";
	for (int m = 0; m <= nmax; m++)
	{
		std::cout << P[nmax+1][m] << '\t' << P[nmax+2][m] << std::endl;
	}
	std::cout << "\033[0m";
	*/

	double v_r, v_theta, v_lambda, tmp1, tmp2;
	double mu = EARTH_MU / r;
	double R = EARTH_RADIUS / r;
	int ind;
	double cnm, snm;

	v_r = 0.0;
	v_theta = 0.0;
	v_lambda = 0.0;

	int N = nmax;

	for (int m = 0; m < nmax+1; m++)
	{
		double cosml = P[N+1][m];
		double sinml = P[N+2][m];
		for (int n = m; n < nmax+1; n++)
		{
			tmp1 = std::pow(R, n);
			tmp2 = tmp1 * (double)(-1.0 * n - 1.0) / r;
			ind = Order2[n][m];
			cnm = Cnm1[ind];
			snm = Snm1[ind];
			
			v_r += mu * (cnm * cosml + snm * sinml) * tmp2 * P[n][m];
			v_theta += mu * (cnm * cosml + snm * sinml) * tmp1 * P1[n][m];
			v_lambda += mu * m * (snm * cosml - cnm * sinml) * tmp1 * P[n][m];
		}
	}

	//std::cout << "\033[32m" << v_r << '\t' << v_theta << '\t' << v_lambda << "\033[0m" << std::endl;

	Result[0] = u * cl * v_r + t * cl / r * v_theta - sl / r / u * v_lambda;
	Result[1] = u * sl * v_r + t * sl / r * v_theta + cl / r / u * v_lambda;
	Result[2] = t * v_r - u / r * v_theta;

	for (int i = 0; i < nmax+3; i++) delete[] P[i];
	delete[] P;

	for (int i = 0; i < nmax+1; i++) delete[] P1[i];
	delete[] P1;

	return;
}