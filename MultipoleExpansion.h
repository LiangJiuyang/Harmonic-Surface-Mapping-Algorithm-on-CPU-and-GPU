#include"InitialSet.h"

//set R^{wavy-line}mn    Now small p=4
int p=15;     


void CalculateMultipleExpansion(double *Q, double x, double y, double z)//Q is a p*p vector.  Set multi-pole expansion coefficient
{
	//double Q[p*p];
	Q[0] = 1.0;
	//QNormlization[0] = Q[0];
	Q[1] = y / 2;
	//QNormlization[1] = Q[1] * sqrt(fac(2)*fac(0));
	Q[2] =(-z);
	//QNormlization[2] = Q[2];
	Q[3] =(-x / 2);
	//QNormlization[3] = Q[3] * sqrt(fac(0)*fac(2));
	//Q[4] = sqrt(fac(4 + 0.0)*fac(0 + 0.0))*(y*Q[2*2 - 1] - x*Q[2*2 - 2 * 2 + 1]) / (2 * 2 + 0.00);
	//Q[4] = sqrt(fac(4)*fac(0))*(y*sqrt(fac(0)*fac(2))*(-x / 2) - x*sqrt(fac(2)*fac(0))*y / 2) / 4;
	int t = 4;
	int m, n;
	//int mm = 4;
	for (int i = 2; i < p; i++)
	{
		while (t<(i+1)*(i+1))
		{
			m = t - i - i*i;
			n = i;
			if (m==-n)
			{
				//Q[t] = (y*Q[n*n - 1] - x*Q[n*n - 2 * n + 1]) / (2 * n + 0.00);//no normalization
				Q[t] =(y*Q[n*n - 1] - x*Q[n*n - 2 * n + 1]) / (2 * n + 0.00);
				//mm++;
			}
			else if (m==n)
			{
				//Q[t] = -(x*Q[n*n - 1] + y*Q[n*n - 2 * n + 1]) / (2 * n + 0.00);//no normlization
				Q[t] = (-(x*Q[n*n - 1] + y*Q[n*n - 2 * n + 1]) / (2 * n + 0.00));
				//mm++;
			}
			else if (n - abs(m) == 1)
			{
				//Q[t] = (-z)*Q[n*n-n+m];//no normlization
				Q[t] = (-z)*Q[n*n-n+m];
				//mm++;
			}
			else if ((n-abs(m))>1)
			{
				//Q[t] = -((2 * n - 1.0)*z*Q[n*n - n + m] + (x*x + y*y + z*z)*Q[n*n - 3 * n + m + 2]) / ((n - abs(m))*(n + abs(m)));//no normlization
				Q[t] = (-((2 * n - 1.0)*z*Q[n*n-n+m] + (x*x + y*y + z*z)*Q[n*n-3*n+m+2]) / ((n - abs(m)+0.0)*(n + abs(m)+0.0)));
				//mm++;
			}
			//QNormlization[t] = Q[t] * sqrt((fac(n - m))*fac(n + m));
			//cout << "Q[" << t << "]=" << Q[t] << endl;
			t++;
		}
	}
	//cout << "mm=" << mm << endl;
	//cout << "t=" << t << endl;
	t = 0;//normlization     Please do not normlize after every step!! That's wrong!
	
	for (int i = 0; i < p; i++)
	{
		while (t < (i + 1)*(i + 1))
		{
			m = t - i - i*i;
			n = i;
			Q[t] = Q[t] * sqrt((fac(n - m))*fac(n + m));
			t++;
		}
	}
	
}
