#include "public.h"
#include<fstream>
#include <iostream>

#include <string>

#include < iomanip >//保留小数点后几位的

#include<cmath>

#include <Eigen/Dense>
#include<Eigen/Core>

const double EARTH_RADIUS = 6371.004;//地球半径km
const double PI = 3.1415926;


using namespace std;






double rad(double d)
{
	return d *  PI/ 180.0;
}


double Get_Distance(double L1, double B1, double L2, double B2)
{
	double radLat1 = rad(B1);
	double radLat2 = rad(B2);
	double a = radLat1 - radLat2;
	double b = rad(L1) - rad(L2);


	double s = 2 * asin(sqrt(pow(sin(a / 2), 2) +
		cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2)));
	s = s * EARTH_RADIUS;
	
	return s;

}




double Median(double p[], int n,bool isAbs)
{


	if (isAbs)
	{
		for (int i = 0; i < n; i++) p[i] = abs(p[i]);
	
	}


	int k = n / 2;
	while (k>0)
	{
		for (int j = k; j < n; j++)
		{
			double t = p[j];
			int i = j - k;
		
			while ((i>=0)&&(p[i]>t))
			{
				p[i + k] = p[i];
				i = i - k;
			}
			p[i + k] = t;
		}
		k = k / 2;
	}
	double mean = (n % 2 == 1) ? p[n / 2] : (p[n / 2] + p[n / 2 - 1]) / 2.0;
	
	return mean;

}
