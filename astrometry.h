#include<./libs/sofa/sofa.h>
#include<./libs/sofa/sofam.h>
#include<string>
#include<array>

int order = 7;
int _mjd = 0;
double mean_earth_rotation_rate = 7.2921158553*1e-5;
double duts[7][2];
double xps[7][2];
double yps[7][2];

double rc2t[3][3];
double rc2t1[3][3];
double rt2c[3][3];

int year = 2022;
int month = 2;
int day = 24;

void EOP()
{
	std::cout << "I'm about to collect EOPs" << std::endl;

	int year_r, month_r, day_r, mjd_r, dat;
	double x_r, y_r, dut_r;
	double lod, dpsi, deps, dx, dy, djm0, djm;
	std::string eopfilename = "../Files/eop_new.txt";
	std::ifstream eop(eopfilename);
	if (!eop.is_open())
	{
		std::cout << "No such EOP file" << std::endl;
	}
	else
	{
		iauCal2jd(year, month, day, &djm0, &djm);
		
		for (int i = 0; i < order; i++)
		{
			eop >> year_r >> month_r >> day_r >> mjd_r >> x_r >> y_r >> dut_r >> lod >> dpsi >> deps >> dx >> dy >> dat;
			duts[i][0] = (double)i * 86400.0;
			xps[i][0] = (double)i * 86400.0;
			yps[i][0] = (double)i * 86400.0;

			duts[i][1] = dut_r;
			xps[i][1] = x_r;
			yps[i][1] = y_r;
		}

		int stop_sign = 0;
		bool found = false;
		while (eop >> year_r >> month_r >> day_r >> mjd_r >> x_r >> y_r >> dut_r >> lod >> dpsi >> deps >> dx >> dy >> dat)
		{
			for (int i = 0; i < order - 1; i++)
			{
				duts[i][1] = duts[i+1][1];
				xps[i][1] = xps[i+1][1];
				yps[i][1] = yps[i+1][1];
			}
			
			duts[order-1][1] = dut_r;
			xps[order-1][1] = x_r;
			yps[order-1][1] = y_r;

			if (mjd_r == (int)djm - order/2)
			{
				found = true;
			}

			if (found)
			{
				//duts[i][0] = (double)i * 86400.0;
				//xps[i][0] = (double)i * 86400.0;
				//yps[i][0] = (double)i * 86400.0;

				//duts[i][1] = dut_r;
				//xps[i][1] = x_r;
				//yps[i][1] = y_r;

				//if (i == 3) _mjd = mjd_r;

				stop_sign++;
			}

			if (stop_sign == order) break;
		}
		eop.close();
	}
}

double lagrange_interpol(double m[][2], double x)
{
	double result = 0.0;
	for (int i = 0; i < order; i++)
	{
		double polinomial = 1.0;
		for (int j = 0; j < order; j++)
		{
			if (j != i)
			{
				polinomial *= (x - m[j][0]) / (m[i][0] - m[j][0]);
			}
		}
		result += m[i][1] * polinomial;
	}
	return result;
}

void rotationMatrices(double seconds)
{
	int j;
	double utc1, utc2, ut11, ut12, tai1, tai2, tt1, tt2;
	double dut, xp, yp;

	double x = (double)(order/2)*86400.0 + seconds; // x is in the middle of the interpolation range
	
	/**/
	dut = lagrange_interpol(duts, x);
	xp = lagrange_interpol(xps, x);
	yp = lagrange_interpol(yps, x);

	//dut = duts[order/2][1];
	//xp = xps[order/2][1];
	//yp = yps[order/2][1];

	xp = xp * M_PI / 648000.0;
	yp = yp * M_PI / 648000.0;

	int hours, minutes;
	hours = (int)(seconds / 3600);
	minutes = (int)((seconds - (double)hours * 3600.0) / 60);

	/*UTC into internal format*/
	j = iauDtf2d("UTC", year, month, day, hours, minutes, seconds - (double)hours * 3600. - (double)minutes * 60.0, &utc1, &utc2);

	/*UTC -> UT1*/
	j = iauUtcut1(utc1, utc2, dut, &ut11, &ut12);

	/*UTC -> TAI*/
	j = iauUtctai(utc1, utc2, &tai1, &tai2);

	/*TAI -> TT*/
	j = iauTaitt(tai1, tai2, &tt1, &tt2);

	iauC2t00b(tt1, tt2, ut11, ut12, xp, yp, rc2t);
	/*
	double sp = iauSp00(tt1, tt2);
	double rpom[3][3];
	iauPom00(xp, yp, sp, rpom);

	double era = iauEra00(ut11, ut12);
	double r[3][3];

	double pnm[3][3];
	iauPnm06a(tt1, tt2, pnm);

	double rerapnm[3][3];
	iauCr(pnm, r);
	iauRz(era, r);
	iauRxr(rpom, r, rc2t_double);
	*/

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			rt2c[i][j] = rc2t[j][i];
		}
	}

	// ---- finding a derivative ---------- //

	double sp = iauSp00(tt1, tt2);
	double rpom[3][3];
	iauPom00(xp, yp, sp, rpom);

	// ��������������� �������
	double m[3][3], thnm[3][3], dthnm[3][3];

	// ���������� Rpm
	iauTr(rpom, rpom);
	iauRxr(rpom, rc2t, thnm);

	// ������ �����������
	m[0][0] = 0.0; m[0][1] = mean_earth_rotation_rate; m[0][2] = 0.0;
	m[1][0] = -1.0 * mean_earth_rotation_rate; m[1][1] = 0.0; m[1][2] = 0.0;
	m[2][0] = 0.0; m[2][1] = 0.0; m[2][2] = 0.0;

	iauRxr(m, thnm, dthnm);

	// ���������� Rpm
	iauTr(rpom, rpom);
	iauRxr(rpom, dthnm, rc2t1);
}