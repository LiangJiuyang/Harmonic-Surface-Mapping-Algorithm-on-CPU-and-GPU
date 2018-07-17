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


/*
double Rm0n0(double x, double y, double z) { return 1; }

double Rm1n1(double x, double y, double z) { return -x/2; }
double Rm_1n1(double x, double y, double z) { return y/2; }
double Rm0n1(double x, double y, double z) { return -z*Rm0n0(x,y,z); }

double Rm2n2(double x, double y, double z) { return -(x*Rm1n1(x,y,z)+y*Rm_1n1(x,y,z))/4; }
double Rm_2n2(double x, double y, double z) { return (y*Rm1n1(x,y,z)-x*Rm_1n1(x,y,z))/4; }
double Rm1n2(double x, double y, double z) { return -z*Rm1n1(x,y,z); }
double Rm_1n2(double x, double y, double z) { return -z*Rm_1n1(x,y,z); }
double Rm0n2(double x, double y, double z) { return -((2*2-1)*z*Rm0n1(x,y,z)+(x*x+y*y+z*z)*Rm0n0(x,y,z))/((2-abs(0))*(2+abs(0))); }

double Rm3n3(double x, double y, double z) { return -(x*Rm2n2(x, y, z) + y*Rm_2n2(x, y, z)) / 6; }
double Rm_3n3(double x, double y, double z) { return (y*Rm2n2(x, y, z) - x*Rm_2n2(x, y, z)) / 6; }
double Rm2n3(double x, double y, double z) { return -z*Rm2n2(x, y, z); }
double Rm_2n3(double x, double y, double z) { return -z*Rm_2n2(x, y, z);}
double Rm1n3(double x, double y, double z) { return -((2*3-1)*z*Rm1n2(x,y,z)+(x*x+y*y+z*z)*Rm1n1(x,y,z))/((3-abs(1))*(3+abs(1))); }
double Rm_1n3(double x, double y, double z) { return -((2 * 3 - 1)*z*Rm_1n2(x, y, z) + (x*x + y*y + z*z)*Rm_1n1(x, y, z)) / ((3 - abs(-1))*(3 + abs(-1))); }
double Rm0n3(double x, double y, double z) { return -((2 * 3 - 1)*z*Rm0n2(x, y, z) + (x*x + y*y + z*z)*Rm0n1(x, y, z)) / ((3 - abs(0))*(3 + abs(0))); }

double Rm4n4(double x, double y, double z) { return -(x*Rm3n3(x, y, z) + y*Rm_3n3(x, y, z)) / 8; }
double Rm_4n4(double x, double y, double z) { return (y*Rm3n3(x, y, z) - x*Rm_3n3(x, y, z)) / 8; }
double Rm3n4(double x, double y, double z) { return -z*Rm3n3(x, y, z); }
double Rm_3n4(double x, double y, double z) { return -z*Rm_3n3(x, y, z); }
double Rm2n4(double x, double y, double z) { return -((2 * 4 - 1)*z*Rm2n3(x, y, z) + (x*x + y*y + z*z)*Rm2n2(x, y, z)) / ((4 - abs(2))*(4 + abs(2))); }
double Rm_2n4(double x, double y, double z) { return -((2 * 4 - 1)*z*Rm_2n3(x, y, z) + (x*x + y*y + z*z)*Rm_2n2(x, y, z)) / ((4 - abs(-2))*(4 + abs(-2))); }
double Rm1n4(double x, double y, double z) { return -((2 * 4 - 1)*z*Rm1n3(x, y, z) + (x*x + y*y + z*z)*Rm1n2(x, y, z)) / ((4 - abs(1))*(4 + abs(1))); }
double Rm_1n4(double x, double y, double z) { return -((2 * 4 - 1)*z*Rm_1n3(x, y, z) + (x*x + y*y + z*z)*Rm_1n2(x, y, z)) / ((4 - abs(-1))*(4 + abs(-1))); }
double Rm0n4(double x, double y, double z) { return -((2 * 4 - 1)*z*Rm0n3(x, y, z) + (x*x + y*y + z*z)*Rm0n2(x, y, z)) / ((4 - abs(0))*(4 + abs(0))); }

double Rm5n5(double x, double y, double z) { return -(x*Rm4n4(x, y, z) + y*Rm_4n4(x, y, z)) / 10; }
double Rm_5n5(double x, double y, double z) { return (y*Rm4n4(x, y, z) - x*Rm_4n4(x, y, z)) / 10; }
double Rm4n5(double x, double y, double z) { return -z*Rm4n4(x, y, z); }
double Rm_4n5(double x, double y, double z) { return -z*Rm_4n4(x, y, z); }
double Rm3n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm3n4(x, y, z) + (x*x + y*y + z*z)*Rm3n3(x, y, z)) / ((5 - abs(3))*(5 + abs(3))); }
double Rm_3n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm_3n4(x, y, z) + (x*x + y*y + z*z)*Rm_3n3(x, y, z)) / ((5 - abs(-3))*(5 + abs(-3))); }
double Rm2n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm2n4(x, y, z) + (x*x + y*y + z*z)*Rm2n3(x, y, z)) / ((5 - abs(2))*(5 + abs(2))); }
double Rm_2n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm_2n4(x, y, z) + (x*x + y*y + z*z)*Rm_2n3(x, y, z)) / ((5 - abs(-2))*(5 + abs(-2))); }
double Rm1n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm1n4(x, y, z) + (x*x + y*y + z*z)*Rm1n3(x, y, z)) / ((5 - abs(1))*(5 + abs(1))); }
double Rm_1n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm_1n4(x, y, z) + (x*x + y*y + z*z)*Rm_1n3(x, y, z)) / ((5 - abs(-1))*(5 + abs(-1))); }
double Rm0n5(double x, double y, double z) { return -((2 * 5 - 1)*z*Rm0n4(x, y, z) + (x*x + y*y + z*z)*Rm0n3(x, y, z)) / ((5 - abs(0))*(5 + abs(0))); }

double Rm6n6(double x, double y, double z) { return -(x*Rm5n5(x, y, z) + y*Rm_5n5(x, y, z)) / 12; }
double Rm_6n6(double x, double y, double z) { return (y*Rm5n5(x, y, z) - x*Rm_5n5(x, y, z)) / 12; }
double Rm5n6(double x, double y, double z) { return -z*Rm5n5(x, y, z); }
double Rm_5n6(double x, double y, double z) { return -z*Rm_5n5(x, y, z); }
double Rm4n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm4n5(x, y, z) + (x*x + y*y + z*z)*Rm4n4(x, y, z)) / ((6 - abs(4))*(6 + abs(4))); }
double Rm_4n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm_4n5(x, y, z) + (x*x + y*y + z*z)*Rm_4n4(x, y, z)) / ((6 - abs(-4))*(6 + abs(-4))); }
double Rm3n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm3n5(x, y, z) + (x*x + y*y + z*z)*Rm3n4(x, y, z)) / ((6 - abs(3))*(6 + abs(3))); }
double Rm_3n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm_3n5(x, y, z) + (x*x + y*y + z*z)*Rm_3n4(x, y, z)) / ((6 - abs(-3))*(6 + abs(-3))); }
double Rm2n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm2n5(x, y, z) + (x*x + y*y + z*z)*Rm2n4(x, y, z)) / ((6 - abs(2))*(6 + abs(2))); }
double Rm_2n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm_2n5(x, y, z) + (x*x + y*y + z*z)*Rm_2n4(x, y, z)) / ((6 - abs(-2))*(6 + abs(-2))); }
double Rm1n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm1n5(x, y, z) + (x*x + y*y + z*z)*Rm1n4(x, y, z)) / ((6 - abs(1))*(6 + abs(1))); }
double Rm_1n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm_1n5(x, y, z) + (x*x + y*y + z*z)*Rm_1n4(x, y, z)) / ((6 - abs(-1))*(6 + abs(-1))); }
double Rm0n6(double x, double y, double z) { return -((2 * 6 - 1)*z*Rm0n5(x, y, z) + (x*x + y*y + z*z)*Rm0n4(x, y, z)) / ((6 - abs(0))*(6 + abs(0))); }

double Rm7n7(double x, double y, double z) { return -(x*Rm6n6(x, y, z) + y*Rm_6n6(x, y, z)) / 14; }
double Rm_7n7(double x, double y, double z) { return (y*Rm6n6(x, y, z) - x*Rm_6n6(x, y, z)) / 14; }
double Rm6n7(double x, double y, double z) { return -z*Rm6n6(x, y, z); }
double Rm_6n7(double x, double y, double z) { return -z*Rm_6n6(x, y, z); }
double Rm5n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm5n6(x, y, z) + (x*x + y*y + z*z)*Rm5n5(x, y, z)) / ((7 - abs(5))*(7 + abs(5))); }
double Rm_5n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm_5n6(x, y, z) + (x*x + y*y + z*z)*Rm_5n5(x, y, z)) / ((7 - abs(-5))*(7 + abs(-5))); }
double Rm4n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm4n6(x, y, z) + (x*x + y*y + z*z)*Rm4n5(x, y, z)) / ((7 - abs(4))*(7 + abs(4))); }
double Rm_4n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm_4n6(x, y, z) + (x*x + y*y + z*z)*Rm_4n5(x, y, z)) / ((7 - abs(-4))*(7 + abs(-4))); }
double Rm3n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm3n6(x, y, z) + (x*x + y*y + z*z)*Rm3n5(x, y, z)) / ((7 - abs(3))*(7 + abs(3))); }
double Rm_3n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm_3n6(x, y, z) + (x*x + y*y + z*z)*Rm_3n5(x, y, z)) / ((7 - abs(-3))*(7 + abs(-3))); }
double Rm2n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm2n6(x, y, z) + (x*x + y*y + z*z)*Rm2n5(x, y, z)) / ((7 - abs(2))*(7 + abs(2))); }
double Rm_2n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm_2n6(x, y, z) + (x*x + y*y + z*z)*Rm_2n5(x, y, z)) / ((7 - abs(-2))*(7 + abs(-2))); }
double Rm1n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm1n6(x, y, z) + (x*x + y*y + z*z)*Rm1n5(x, y, z)) / ((7 - abs(1))*(7 + abs(1))); }
double Rm_1n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm_1n6(x, y, z) + (x*x + y*y + z*z)*Rm_1n5(x, y, z)) / ((7 - abs(-1))*(7 + abs(-1))); }
double Rm0n7(double x, double y, double z) { return -((2 * 7 - 1)*z*Rm0n6(x, y, z) + (x*x + y*y + z*z)*Rm0n5(x, y, z)) / ((7 - abs(0))*(7 + abs(0))); }

double Rm8n8(double x, double y, double z) { return -(x*Rm7n7(x, y, z) + y*Rm_7n7(x, y, z)) / 16; }
double Rm_8n8(double x, double y, double z) { return (y*Rm7n7(x, y, z) - x*Rm_7n7(x, y, z)) / 16; }
double Rm7n8(double x, double y, double z) { return -z*Rm7n7(x, y, z); }
double Rm_7n8(double x, double y, double z) { return -z*Rm_7n7(x, y, z); }
double Rm6n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm6n7(x, y, z) + (x*x + y*y + z*z)*Rm6n6(x, y, z)) / ((8 - abs(6))*(8 + abs(6))); }
double Rm_6n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_6n7(x, y, z) + (x*x + y*y + z*z)*Rm_6n6(x, y, z)) / ((8 - abs(-6))*(8 + abs(-6))); }
double Rm5n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm5n7(x, y, z) + (x*x + y*y + z*z)*Rm5n6(x, y, z)) / ((8 - abs(5))*(8 + abs(5))); }
double Rm_5n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_5n7(x, y, z) + (x*x + y*y + z*z)*Rm_5n6(x, y, z)) / ((8 - abs(-5))*(8 + abs(-5))); }
double Rm4n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm4n7(x, y, z) + (x*x + y*y + z*z)*Rm4n6(x, y, z)) / ((8 - abs(4))*(8 + abs(4))); }
double Rm_4n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_4n7(x, y, z) + (x*x + y*y + z*z)*Rm_4n6(x, y, z)) / ((8 - abs(-4))*(8 + abs(-4))); }
double Rm3n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm3n7(x, y, z) + (x*x + y*y + z*z)*Rm3n6(x, y, z)) / ((8 - abs(3))*(8 + abs(3))); }
double Rm_3n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_3n7(x, y, z) + (x*x + y*y + z*z)*Rm_3n6(x, y, z)) / ((8 - abs(-3))*(8 + abs(-3))); }
double Rm2n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm2n7(x, y, z) + (x*x + y*y + z*z)*Rm2n6(x, y, z)) / ((8 - abs(2))*(8 + abs(2))); }
double Rm_2n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_2n7(x, y, z) + (x*x + y*y + z*z)*Rm_2n6(x, y, z)) / ((8 - abs(-2))*(8 + abs(-2))); }
double Rm1n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm1n7(x, y, z) + (x*x + y*y + z*z)*Rm1n6(x, y, z)) / ((8 - abs(1))*(8 + abs(1))); }
double Rm_1n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm_1n7(x, y, z) + (x*x + y*y + z*z)*Rm_1n6(x, y, z)) / ((8 - abs(-1))*(8 + abs(-1))); }
double Rm0n8(double x, double y, double z) { return -((2 * 8 - 1)*z*Rm0n7(x, y, z) + (x*x + y*y + z*z)*Rm0n6(x, y, z)) / ((8 - abs(0))*(8 + abs(0))); }

double Rm9n9(double x, double y, double z) { return -(x*Rm8n8(x, y, z) + y*Rm_8n8(x, y, z)) / 18; }
double Rm_9n9(double x, double y, double z) { return (y*Rm8n8(x, y, z) - x*Rm_8n8(x, y, z)) / 18; }
double Rm8n9(double x, double y, double z) { return -z*Rm8n8(x, y, z); }
double Rm_8n9(double x, double y, double z) { return -z*Rm_8n8(x, y, z); }
double Rm7n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm7n8(x, y, z) + (x*x + y*y + z*z)*Rm7n7(x, y, z)) / ((9 - abs(7))*(9 + abs(7))); }
double Rm_7n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_7n8(x, y, z) + (x*x + y*y + z*z)*Rm_7n7(x, y, z)) / ((9 - abs(-7))*(9 + abs(-7))); }
double Rm6n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm6n8(x, y, z) + (x*x + y*y + z*z)*Rm6n7(x, y, z)) / ((9 - abs(6))*(9 + abs(6))); }
double Rm_6n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_6n8(x, y, z) + (x*x + y*y + z*z)*Rm_6n7(x, y, z)) / ((9 - abs(-6))*(9 + abs(-6))); }
double Rm5n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm5n8(x, y, z) + (x*x + y*y + z*z)*Rm5n7(x, y, z)) / ((9 - abs(5))*(9 + abs(5))); }
double Rm_5n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_5n8(x, y, z) + (x*x + y*y + z*z)*Rm_5n7(x, y, z)) / ((9 - abs(-5))*(9 + abs(-5))); }
double Rm4n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm4n8(x, y, z) + (x*x + y*y + z*z)*Rm4n7(x, y, z)) / ((9 - abs(4))*(9 + abs(4))); }
double Rm_4n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_4n8(x, y, z) + (x*x + y*y + z*z)*Rm_4n7(x, y, z)) / ((9 - abs(-4))*(9 + abs(-4))); }
double Rm3n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm3n8(x, y, z) + (x*x + y*y + z*z)*Rm3n7(x, y, z)) / ((9 - abs(3))*(9 + abs(3))); }
double Rm_3n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_3n8(x, y, z) + (x*x + y*y + z*z)*Rm_3n7(x, y, z)) / ((9 - abs(-3))*(9 + abs(-3))); }
double Rm2n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm2n8(x, y, z) + (x*x + y*y + z*z)*Rm2n7(x, y, z)) / ((9 - abs(2))*(9 + abs(2))); }
double Rm_2n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_2n8(x, y, z) + (x*x + y*y + z*z)*Rm_2n7(x, y, z)) / ((9 - abs(-2))*(9 + abs(-2))); }
double Rm1n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm1n8(x, y, z) + (x*x + y*y + z*z)*Rm1n7(x, y, z)) / ((9 - abs(1))*(9 + abs(1))); }
double Rm_1n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm_1n8(x, y, z) + (x*x + y*y + z*z)*Rm_1n7(x, y, z)) / ((9 - abs(-1))*(9 + abs(-1))); }
double Rm0n9(double x, double y, double z) { return -((2 * 9 - 1)*z*Rm0n8(x, y, z) + (x*x + y*y + z*z)*Rm0n7(x, y, z)) / ((9 - abs(0))*(9 + abs(0))); }

double Rm10n10(double x, double y, double z) { return -(x*Rm9n9(x, y, z) + y*Rm_9n9(x, y, z)) / 20; }
double Rm_10n10(double x, double y, double z) { return (y*Rm9n9(x, y, z) - x*Rm_9n9(x, y, z)) / 20; }
double Rm9n10(double x, double y, double z) { return -z*Rm9n9(x, y, z); }
double Rm_9n10(double x, double y, double z) { return -z*Rm_9n9(x, y, z); }
double Rm8n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm8n9(x, y, z) + (x*x + y*y + z*z)*Rm8n8(x, y, z)) / ((10 - abs(8))*(10 + abs(8))); }
double Rm_8n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_8n9(x, y, z) + (x*x + y*y + z*z)*Rm_8n8(x, y, z)) / ((10 - abs(-8))*(10 + abs(-8))); }
double Rm7n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm7n9(x, y, z) + (x*x + y*y + z*z)*Rm7n8(x, y, z)) / ((10 - abs(7))*(10 + abs(7))); }
double Rm_7n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_7n9(x, y, z) + (x*x + y*y + z*z)*Rm_7n8(x, y, z)) / ((10 - abs(-7))*(10 + abs(-7))); }
double Rm6n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm6n9(x, y, z) + (x*x + y*y + z*z)*Rm6n8(x, y, z)) / ((10 - abs(6))*(10 + abs(6))); }
double Rm_6n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_6n9(x, y, z) + (x*x + y*y + z*z)*Rm_6n8(x, y, z)) / ((10 - abs(-6))*(10 + abs(-6))); }
double Rm5n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm5n9(x, y, z) + (x*x + y*y + z*z)*Rm5n8(x, y, z)) / ((10 - abs(5))*(10 + abs(5))); }
double Rm_5n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_5n9(x, y, z) + (x*x + y*y + z*z)*Rm_5n8(x, y, z)) / ((10 - abs(-5))*(10 + abs(-5))); }
double Rm4n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm4n9(x, y, z) + (x*x + y*y + z*z)*Rm4n8(x, y, z)) / ((10 - abs(4))*(10 + abs(4))); }
double Rm_4n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_4n9(x, y, z) + (x*x + y*y + z*z)*Rm_4n8(x, y, z)) / ((10 - abs(-4))*(10 + abs(-4))); }
double Rm3n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm3n9(x, y, z) + (x*x + y*y + z*z)*Rm3n8(x, y, z)) / ((10 - abs(3))*(10 + abs(3))); }
double Rm_3n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_3n9(x, y, z) + (x*x + y*y + z*z)*Rm_3n8(x, y, z)) / ((10 - abs(-3))*(10 + abs(-3))); }
double Rm2n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm2n9(x, y, z) + (x*x + y*y + z*z)*Rm2n8(x, y, z)) / ((10 - abs(2))*(10 + abs(2))); }
double Rm_2n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_2n9(x, y, z) + (x*x + y*y + z*z)*Rm_2n8(x, y, z)) / ((10 - abs(-2))*(10 + abs(-2))); }
double Rm1n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm1n9(x, y, z) + (x*x + y*y + z*z)*Rm1n8(x, y, z)) / ((10 - abs(1))*(10 + abs(1))); }
double Rm_1n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm_1n9(x, y, z) + (x*x + y*y + z*z)*Rm_1n8(x, y, z)) / ((10 - abs(-1))*(10 + abs(-1))); }
double Rm0n10(double x, double y, double z) { return -((2 * 10 - 1)*z*Rm0n9(x, y, z) + (x*x + y*y + z*z)*Rm0n8(x, y, z)) / ((10 - abs(0))*(10 + abs(0))); }

double Rm11n11(double x, double y, double z) { return -(x*Rm10n10(x, y, z) + y*Rm_10n10(x, y, z)) / 22; }
double Rm_11n11(double x, double y, double z) { return (y*Rm10n10(x, y, z) - x*Rm_10n10(x, y, z)) / 22; }
double Rm10n11(double x, double y, double z) { return -z*Rm10n10(x, y, z); }
double Rm_10n11(double x, double y, double z) { return -z*Rm_10n10(x, y, z); }
double Rm9n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm9n10(x, y, z) + (x*x + y*y + z*z)*Rm9n9(x, y, z)) / ((11 - abs(9))*(11 + abs(9))); }
double Rm_9n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_9n10(x, y, z) + (x*x + y*y + z*z)*Rm_9n9(x, y, z)) / ((11 - abs(-9))*(11 + abs(-9))); }
double Rm8n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm8n10(x, y, z) + (x*x + y*y + z*z)*Rm8n9(x, y, z)) / ((11 - abs(8))*(11 + abs(8))); }
double Rm_8n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_8n10(x, y, z) + (x*x + y*y + z*z)*Rm_8n9(x, y, z)) / ((11 - abs(-8))*(11 + abs(-8))); }
double Rm7n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm7n10(x, y, z) + (x*x + y*y + z*z)*Rm7n9(x, y, z)) / ((11 - abs(7))*(11 + abs(7))); }
double Rm_7n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_7n10(x, y, z) + (x*x + y*y + z*z)*Rm_7n9(x, y, z)) / ((11 - abs(-7))*(11 + abs(-7))); }
double Rm6n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm6n10(x, y, z) + (x*x + y*y + z*z)*Rm6n9(x, y, z)) / ((11 - abs(6))*(11 + abs(6))); }
double Rm_6n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_6n10(x, y, z) + (x*x + y*y + z*z)*Rm_6n9(x, y, z)) / ((11 - abs(-6))*(11 + abs(-6))); }
double Rm5n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm5n10(x, y, z) + (x*x + y*y + z*z)*Rm5n9(x, y, z)) / ((11 - abs(5))*(11 + abs(5))); }
double Rm_5n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_5n10(x, y, z) + (x*x + y*y + z*z)*Rm_5n9(x, y, z)) / ((11 - abs(-5))*(11 + abs(-5))); }
double Rm4n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm4n10(x, y, z) + (x*x + y*y + z*z)*Rm4n9(x, y, z)) / ((11 - abs(4))*(11 + abs(4))); }
double Rm_4n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_4n10(x, y, z) + (x*x + y*y + z*z)*Rm_4n9(x, y, z)) / ((11 - abs(-4))*(11 + abs(-4))); }
double Rm3n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm3n10(x, y, z) + (x*x + y*y + z*z)*Rm3n9(x, y, z)) / ((11 - abs(3))*(11 + abs(3))); }
double Rm_3n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_3n10(x, y, z) + (x*x + y*y + z*z)*Rm_3n9(x, y, z)) / ((11 - abs(-3))*(11 + abs(-3))); }
double Rm2n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm2n10(x, y, z) + (x*x + y*y + z*z)*Rm2n9(x, y, z)) / ((11 - abs(2))*(11 + abs(2))); }
double Rm_2n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_2n10(x, y, z) + (x*x + y*y + z*z)*Rm_2n9(x, y, z)) / ((11 - abs(-2))*(11 + abs(-2))); }
double Rm1n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm1n10(x, y, z) + (x*x + y*y + z*z)*Rm1n9(x, y, z)) / ((11 - abs(1))*(11 + abs(1))); }
double Rm_1n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm_1n10(x, y, z) + (x*x + y*y + z*z)*Rm_1n9(x, y, z)) / ((11 - abs(-1))*(11 + abs(-1))); }
double Rm0n11(double x, double y, double z) { return -((2 * 11 - 1)*z*Rm0n10(x, y, z) + (x*x + y*y + z*z)*Rm0n9(x, y, z)) / ((11 - abs(0))*(11 + abs(0))); }

double Rm12n12(double x, double y, double z) { return -(x*Rm11n11(x, y, z) + y*Rm_11n11(x, y, z)) / 24; }
double Rm_12n12(double x, double y, double z) { return (y*Rm11n11(x, y, z) - x*Rm_11n11(x, y, z)) / 24; }
double Rm11n12(double x, double y, double z) { return -z*Rm11n11(x, y, z); }
double Rm_11n12(double x, double y, double z) { return -z*Rm_11n11(x, y, z); }
double Rm10n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm10n11(x, y, z) + (x*x + y*y + z*z)*Rm10n10(x, y, z)) / ((12 - abs(10))*(12 + abs(10))); }
double Rm_10n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_10n11(x, y, z) + (x*x + y*y + z*z)*Rm_10n10(x, y, z)) / ((12 - abs(-10))*(12 + abs(-10))); }
double Rm9n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm9n11(x, y, z) + (x*x + y*y + z*z)*Rm9n10(x, y, z)) / ((12 - abs(9))*(12 + abs(9))); }
double Rm_9n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_9n11(x, y, z) + (x*x + y*y + z*z)*Rm_9n10(x, y, z)) / ((12 - abs(-9))*(12 + abs(-9))); }
double Rm8n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm8n11(x, y, z) + (x*x + y*y + z*z)*Rm8n10(x, y, z)) / ((12 - abs(8))*(12 + abs(8))); }
double Rm_8n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_8n11(x, y, z) + (x*x + y*y + z*z)*Rm_8n10(x, y, z)) / ((12 - abs(-8))*(12 + abs(-8))); }
double Rm7n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm7n11(x, y, z) + (x*x + y*y + z*z)*Rm7n10(x, y, z)) / ((12 - abs(7))*(12 + abs(7))); }
double Rm_7n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_7n11(x, y, z) + (x*x + y*y + z*z)*Rm_7n10(x, y, z)) / ((12 - abs(-7))*(12 + abs(-7))); }
double Rm6n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm6n11(x, y, z) + (x*x + y*y + z*z)*Rm6n10(x, y, z)) / ((12 - abs(6))*(12 + abs(6))); }
double Rm_6n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_6n11(x, y, z) + (x*x + y*y + z*z)*Rm_6n10(x, y, z)) / ((12 - abs(-6))*(12 + abs(-6))); }
double Rm5n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm5n11(x, y, z) + (x*x + y*y + z*z)*Rm5n10(x, y, z)) / ((12 - abs(5))*(12 + abs(5))); }
double Rm_5n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_5n11(x, y, z) + (x*x + y*y + z*z)*Rm_5n10(x, y, z)) / ((12 - abs(-5))*(12 + abs(-5))); }
double Rm4n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm4n11(x, y, z) + (x*x + y*y + z*z)*Rm4n10(x, y, z)) / ((12 - abs(4))*(12 + abs(4))); }
double Rm_4n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_4n11(x, y, z) + (x*x + y*y + z*z)*Rm_4n10(x, y, z)) / ((12 - abs(-4))*(12 + abs(-4))); }
double Rm3n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm3n11(x, y, z) + (x*x + y*y + z*z)*Rm3n10(x, y, z)) / ((12 - abs(3))*(12 + abs(3))); }
double Rm_3n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_3n11(x, y, z) + (x*x + y*y + z*z)*Rm_3n10(x, y, z)) / ((12 - abs(-3))*(12 + abs(-3))); }
double Rm2n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm2n11(x, y, z) + (x*x + y*y + z*z)*Rm2n10(x, y, z)) / ((12 - abs(2))*(12 + abs(2))); }
double Rm_2n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_2n11(x, y, z) + (x*x + y*y + z*z)*Rm_2n10(x, y, z)) / ((12 - abs(-2))*(12 + abs(-2))); }
double Rm1n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm1n11(x, y, z) + (x*x + y*y + z*z)*Rm1n10(x, y, z)) / ((12 - abs(1))*(12 + abs(1))); }
double Rm_1n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm_1n11(x, y, z) + (x*x + y*y + z*z)*Rm_1n10(x, y, z)) / ((12 - abs(-1))*(12 + abs(-1))); }
double Rm0n12(double x, double y, double z) { return -((2 * 12 - 1)*z*Rm0n11(x, y, z) + (x*x + y*y + z*z)*Rm0n10(x, y, z)) / ((12 - abs(0))*(12 + abs(0))); }

double Rm13n13(double x, double y, double z) { return -(x*Rm12n12(x, y, z) + y*Rm_12n12(x, y, z)) / 26; }
double Rm_13n13(double x, double y, double z) { return (y*Rm12n12(x, y, z) - x*Rm_12n12(x, y, z)) / 26; }
double Rm12n13(double x, double y, double z) { return -z*Rm12n12(x, y, z); }
double Rm_12n13(double x, double y, double z) { return -z*Rm_12n12(x, y, z); }
double Rm11n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm11n12(x, y, z) + (x*x + y*y + z*z)*Rm11n11(x, y, z)) / ((13 - abs(11))*(13 + abs(11))); }
double Rm_11n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_11n12(x, y, z) + (x*x + y*y + z*z)*Rm_11n11(x, y, z)) / ((13 - abs(-11))*(13 + abs(-11))); }
double Rm10n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm10n12(x, y, z) + (x*x + y*y + z*z)*Rm10n11(x, y, z)) / ((13 - abs(10))*(13 + abs(10))); }
double Rm_10n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_10n12(x, y, z) + (x*x + y*y + z*z)*Rm_10n11(x, y, z)) / ((13 - abs(-10))*(13 + abs(-10))); }
double Rm9n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm9n12(x, y, z) + (x*x + y*y + z*z)*Rm9n11(x, y, z)) / ((13 - abs(9))*(13 + abs(9))); }
double Rm_9n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_9n12(x, y, z) + (x*x + y*y + z*z)*Rm_9n11(x, y, z)) / ((13 - abs(-9))*(13 + abs(-9))); }
double Rm8n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm8n12(x, y, z) + (x*x + y*y + z*z)*Rm8n11(x, y, z)) / ((13 - abs(8))*(13 + abs(8))); }
double Rm_8n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_8n12(x, y, z) + (x*x + y*y + z*z)*Rm_8n11(x, y, z)) / ((13 - abs(-8))*(13 + abs(-8))); }
double Rm7n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm7n12(x, y, z) + (x*x + y*y + z*z)*Rm7n11(x, y, z)) / ((13 - abs(7))*(13 + abs(7))); }
double Rm_7n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_7n12(x, y, z) + (x*x + y*y + z*z)*Rm_7n11(x, y, z)) / ((13 - abs(-7))*(13 + abs(-7))); }
double Rm6n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm6n12(x, y, z) + (x*x + y*y + z*z)*Rm6n11(x, y, z)) / ((13 - abs(6))*(13 + abs(6))); }
double Rm_6n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_6n12(x, y, z) + (x*x + y*y + z*z)*Rm_6n11(x, y, z)) / ((13 - abs(-6))*(13 + abs(-6))); }
double Rm5n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm5n12(x, y, z) + (x*x + y*y + z*z)*Rm5n11(x, y, z)) / ((13 - abs(5))*(13 + abs(5))); }
double Rm_5n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_5n12(x, y, z) + (x*x + y*y + z*z)*Rm_5n11(x, y, z)) / ((13 - abs(-5))*(13 + abs(-5))); }
double Rm4n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm4n12(x, y, z) + (x*x + y*y + z*z)*Rm4n11(x, y, z)) / ((13 - abs(4))*(13 + abs(4))); }
double Rm_4n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_4n12(x, y, z) + (x*x + y*y + z*z)*Rm_4n11(x, y, z)) / ((13 - abs(-4))*(13 + abs(-4))); }
double Rm3n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm3n12(x, y, z) + (x*x + y*y + z*z)*Rm3n11(x, y, z)) / ((13 - abs(3))*(13 + abs(3))); }
double Rm_3n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_3n12(x, y, z) + (x*x + y*y + z*z)*Rm_3n11(x, y, z)) / ((13 - abs(-3))*(13 + abs(-3))); }
double Rm2n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm2n12(x, y, z) + (x*x + y*y + z*z)*Rm2n11(x, y, z)) / ((13 - abs(2))*(13 + abs(2))); }
double Rm_2n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_2n12(x, y, z) + (x*x + y*y + z*z)*Rm_2n11(x, y, z)) / ((13 - abs(-2))*(13 + abs(-2))); }
double Rm1n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm1n12(x, y, z) + (x*x + y*y + z*z)*Rm1n11(x, y, z)) / ((13 - abs(1))*(13 + abs(1))); }
double Rm_1n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm_1n12(x, y, z) + (x*x + y*y + z*z)*Rm_1n11(x, y, z)) / ((13 - abs(-1))*(13 + abs(-1))); }
double Rm0n13(double x, double y, double z) { return -((2 * 13 - 1)*z*Rm0n12(x, y, z) + (x*x + y*y + z*z)*Rm0n11(x, y, z)) / ((13 - abs(0))*(13 + abs(0))); }

double Rm14n14(double x, double y, double z) { return -(x*Rm13n13(x, y, z) + y*Rm_13n13(x, y, z)) / 28; }
double Rm_14n14(double x, double y, double z) { return (y*Rm13n13(x, y, z) - x*Rm_13n13(x, y, z)) / 28; }
double Rm13n14(double x, double y, double z) { return -z*Rm13n13(x, y, z); }
double Rm_13n14(double x, double y, double z) { return -z*Rm_13n13(x, y, z); }
double Rm12n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm12n13(x, y, z) + (x*x + y*y + z*z)*Rm12n12(x, y, z)) / ((14 - abs(12))*(14 + abs(12))); }
double Rm_12n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_12n13(x, y, z) + (x*x + y*y + z*z)*Rm_12n12(x, y, z)) / ((14 - abs(-12))*(14 + abs(-12))); }
double Rm11n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm11n13(x, y, z) + (x*x + y*y + z*z)*Rm11n12(x, y, z)) / ((14 - abs(11))*(14 + abs(11))); }
double Rm_11n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_11n13(x, y, z) + (x*x + y*y + z*z)*Rm_11n12(x, y, z)) / ((14 - abs(-11))*(14 + abs(-11))); }
double Rm10n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm10n13(x, y, z) + (x*x + y*y + z*z)*Rm10n12(x, y, z)) / ((14 - abs(10))*(14 + abs(10))); }
double Rm_10n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_10n13(x, y, z) + (x*x + y*y + z*z)*Rm_10n12(x, y, z)) / ((14 - abs(-10))*(14 + abs(-10))); }
double Rm9n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm9n13(x, y, z) + (x*x + y*y + z*z)*Rm9n12(x, y, z)) / ((14 - abs(9))*(14 + abs(9))); }
double Rm_9n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_9n13(x, y, z) + (x*x + y*y + z*z)*Rm_9n12(x, y, z)) / ((14 - abs(-9))*(14 + abs(-9))); }
double Rm8n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm8n13(x, y, z) + (x*x + y*y + z*z)*Rm8n12(x, y, z)) / ((14 - abs(8))*(14 + abs(8))); }
double Rm_8n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_8n13(x, y, z) + (x*x + y*y + z*z)*Rm_8n12(x, y, z)) / ((14 - abs(-8))*(14 + abs(-8))); }
double Rm7n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm7n13(x, y, z) + (x*x + y*y + z*z)*Rm7n12(x, y, z)) / ((14 - abs(7))*(14 + abs(7))); }
double Rm_7n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_7n13(x, y, z) + (x*x + y*y + z*z)*Rm_7n12(x, y, z)) / ((14 - abs(-7))*(14 + abs(-7))); }
double Rm6n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm6n13(x, y, z) + (x*x + y*y + z*z)*Rm6n12(x, y, z)) / ((14 - abs(6))*(14 + abs(6))); }
double Rm_6n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_6n13(x, y, z) + (x*x + y*y + z*z)*Rm_6n12(x, y, z)) / ((14 - abs(-6))*(14 + abs(-6))); }
double Rm5n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm5n13(x, y, z) + (x*x + y*y + z*z)*Rm5n12(x, y, z)) / ((14 - abs(5))*(14 + abs(5))); }
double Rm_5n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_5n13(x, y, z) + (x*x + y*y + z*z)*Rm_5n12(x, y, z)) / ((14 - abs(-5))*(14 + abs(-5))); }
double Rm4n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm4n13(x, y, z) + (x*x + y*y + z*z)*Rm4n12(x, y, z)) / ((14 - abs(4))*(14 + abs(4))); }
double Rm_4n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_4n13(x, y, z) + (x*x + y*y + z*z)*Rm_4n12(x, y, z)) / ((14 - abs(-4))*(14 + abs(-4))); }
double Rm3n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm3n13(x, y, z) + (x*x + y*y + z*z)*Rm3n12(x, y, z)) / ((14 - abs(3))*(14 + abs(3))); }
double Rm_3n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_3n13(x, y, z) + (x*x + y*y + z*z)*Rm_3n12(x, y, z)) / ((14 - abs(-3))*(14 + abs(-3))); }
double Rm2n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm2n13(x, y, z) + (x*x + y*y + z*z)*Rm2n12(x, y, z)) / ((14 - abs(2))*(14 + abs(2))); }
double Rm_2n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_2n13(x, y, z) + (x*x + y*y + z*z)*Rm_2n12(x, y, z)) / ((14 - abs(-2))*(14 + abs(-2))); }
double Rm1n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm1n13(x, y, z) + (x*x + y*y + z*z)*Rm1n12(x, y, z)) / ((14 - abs(1))*(14 + abs(1))); }
double Rm_1n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm_1n13(x, y, z) + (x*x + y*y + z*z)*Rm_1n12(x, y, z)) / ((14 - abs(-1))*(14 + abs(-1))); }
double Rm0n14(double x, double y, double z) { return -((2 * 14 - 1)*z*Rm0n13(x, y, z) + (x*x + y*y + z*z)*Rm0n12(x, y, z)) / ((14 - abs(0))*(14 + abs(0))); }

double Rm15n15(double x, double y, double z) { return -(x*Rm14n14(x, y, z) + y*Rm_14n14(x, y, z)) / 30; }
double Rm_15n15(double x, double y, double z) { return (y*Rm14n14(x, y, z) - x*Rm_14n14(x, y, z)) / 30; }
double Rm14n15(double x, double y, double z) { return -z*Rm14n14(x, y, z); }
double Rm_14n15(double x, double y, double z) { return -z*Rm_14n14(x, y, z); }
double Rm13n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm13n14(x, y, z) + (x*x + y*y + z*z)*Rm13n13(x, y, z)) / ((15 - abs(13))*(15 + abs(13))); }
double Rm_13n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_13n14(x, y, z) + (x*x + y*y + z*z)*Rm_13n13(x, y, z)) / ((15 - abs(-13))*(15 + abs(-13))); }
double Rm12n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm12n14(x, y, z) + (x*x + y*y + z*z)*Rm12n13(x, y, z)) / ((15 - abs(12))*(15 + abs(12))); }
double Rm_12n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_12n14(x, y, z) + (x*x + y*y + z*z)*Rm_12n13(x, y, z)) / ((15 - abs(-12))*(15 + abs(-12))); }
double Rm11n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm11n14(x, y, z) + (x*x + y*y + z*z)*Rm11n13(x, y, z)) / ((15 - abs(11))*(15 + abs(11))); }
double Rm_11n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_11n14(x, y, z) + (x*x + y*y + z*z)*Rm_11n13(x, y, z)) / ((15 - abs(-11))*(15 + abs(-11))); }
double Rm10n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm10n14(x, y, z) + (x*x + y*y + z*z)*Rm10n13(x, y, z)) / ((15 - abs(10))*(15 + abs(10))); }
double Rm_10n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_10n14(x, y, z) + (x*x + y*y + z*z)*Rm_10n13(x, y, z)) / ((15 - abs(-10))*(15 + abs(-10))); }
double Rm9n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm9n14(x, y, z) + (x*x + y*y + z*z)*Rm9n13(x, y, z)) / ((15 - abs(9))*(15 + abs(9))); }
double Rm_9n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_9n14(x, y, z) + (x*x + y*y + z*z)*Rm_9n13(x, y, z)) / ((15 - abs(-9))*(15 + abs(-9))); }
double Rm8n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm8n14(x, y, z) + (x*x + y*y + z*z)*Rm8n13(x, y, z)) / ((15 - abs(8))*(15 + abs(8))); }
double Rm_8n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_8n14(x, y, z) + (x*x + y*y + z*z)*Rm_8n13(x, y, z)) / ((15 - abs(-8))*(15 + abs(-8))); }
double Rm7n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm7n14(x, y, z) + (x*x + y*y + z*z)*Rm7n13(x, y, z)) / ((15 - abs(7))*(15 + abs(7))); }
double Rm_7n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_7n14(x, y, z) + (x*x + y*y + z*z)*Rm_7n13(x, y, z)) / ((15 - abs(-7))*(15 + abs(-7))); }
double Rm6n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm6n14(x, y, z) + (x*x + y*y + z*z)*Rm6n13(x, y, z)) / ((15 - abs(6))*(15 + abs(6))); }
double Rm_6n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_6n14(x, y, z) + (x*x + y*y + z*z)*Rm_6n13(x, y, z)) / ((15 - abs(-6))*(15 + abs(-6))); }
double Rm5n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm5n14(x, y, z) + (x*x + y*y + z*z)*Rm5n13(x, y, z)) / ((15 - abs(5))*(15 + abs(5))); }
double Rm_5n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_5n14(x, y, z) + (x*x + y*y + z*z)*Rm_5n13(x, y, z)) / ((15 - abs(-5))*(15 + abs(-5))); }
double Rm4n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm4n14(x, y, z) + (x*x + y*y + z*z)*Rm4n13(x, y, z)) / ((15 - abs(4))*(15 + abs(4))); }
double Rm_4n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_4n14(x, y, z) + (x*x + y*y + z*z)*Rm_4n13(x, y, z)) / ((15 - abs(-4))*(15 + abs(-4))); }
double Rm3n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm3n14(x, y, z) + (x*x + y*y + z*z)*Rm3n13(x, y, z)) / ((15 - abs(3))*(15 + abs(3))); }
double Rm_3n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_3n14(x, y, z) + (x*x + y*y + z*z)*Rm_3n13(x, y, z)) / ((15 - abs(-3))*(15 + abs(-3))); }
double Rm2n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm2n14(x, y, z) + (x*x + y*y + z*z)*Rm2n13(x, y, z)) / ((15 - abs(2))*(15 + abs(2))); }
double Rm_2n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_2n14(x, y, z) + (x*x + y*y + z*z)*Rm_2n13(x, y, z)) / ((15 - abs(-2))*(15 + abs(-2))); }
double Rm1n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm1n14(x, y, z) + (x*x + y*y + z*z)*Rm1n13(x, y, z)) / ((15 - abs(1))*(15 + abs(1))); }
double Rm_1n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm_1n14(x, y, z) + (x*x + y*y + z*z)*Rm_1n13(x, y, z)) / ((15 - abs(-1))*(15 + abs(-1))); }
double Rm0n15(double x, double y, double z) { return -((2 * 15 - 1)*z*Rm0n14(x, y, z) + (x*x + y*y + z*z)*Rm0n13(x, y, z)) / ((15 - abs(0))*(15 + abs(0))); }

double Rm16n16(double x, double y, double z) { return -(x*Rm15n15(x, y, z) + y*Rm_15n15(x, y, z)) / 32; }
double Rm_16n16(double x, double y, double z) { return (y*Rm15n15(x, y, z) - x*Rm_15n15(x, y, z)) / 32; }
double Rm15n16(double x, double y, double z) { return -z*Rm15n15(x, y, z); }
double Rm_15n16(double x, double y, double z) { return -z*Rm_15n15(x, y, z); }
double Rm14n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm14n15(x, y, z) + (x*x + y*y + z*z)*Rm14n14(x, y, z)) / ((16 - abs(14))*(16 + abs(14))); }
double Rm_14n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_14n15(x, y, z) + (x*x + y*y + z*z)*Rm_14n14(x, y, z)) / ((16 - abs(-14))*(16 + abs(-14))); }
double Rm13n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm13n15(x, y, z) + (x*x + y*y + z*z)*Rm13n14(x, y, z)) / ((16 - abs(13))*(16 + abs(13))); }
double Rm_13n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_13n15(x, y, z) + (x*x + y*y + z*z)*Rm_13n14(x, y, z)) / ((16 - abs(-13))*(16 + abs(-13))); }
double Rm12n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm12n15(x, y, z) + (x*x + y*y + z*z)*Rm12n14(x, y, z)) / ((16 - abs(12))*(16 + abs(12))); }
double Rm_12n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_12n15(x, y, z) + (x*x + y*y + z*z)*Rm_12n14(x, y, z)) / ((16 - abs(-12))*(16 + abs(-12))); }
double Rm11n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm11n15(x, y, z) + (x*x + y*y + z*z)*Rm11n14(x, y, z)) / ((16 - abs(11))*(16 + abs(11))); }
double Rm_11n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_11n15(x, y, z) + (x*x + y*y + z*z)*Rm_11n14(x, y, z)) / ((16 - abs(-11))*(16 + abs(-11))); }
double Rm10n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm10n15(x, y, z) + (x*x + y*y + z*z)*Rm10n14(x, y, z)) / ((16 - abs(10))*(16 + abs(10))); }
double Rm_10n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_10n15(x, y, z) + (x*x + y*y + z*z)*Rm_10n14(x, y, z)) / ((16 - abs(-10))*(16 + abs(-10))); }
double Rm9n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm9n15(x, y, z) + (x*x + y*y + z*z)*Rm9n14(x, y, z)) / ((16 - abs(9))*(16 + abs(9))); }
double Rm_9n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_9n15(x, y, z) + (x*x + y*y + z*z)*Rm_9n14(x, y, z)) / ((16 - abs(-9))*(16 + abs(-9))); }
double Rm8n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm8n15(x, y, z) + (x*x + y*y + z*z)*Rm8n14(x, y, z)) / ((16 - abs(8))*(16 + abs(8))); }
double Rm_8n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_8n15(x, y, z) + (x*x + y*y + z*z)*Rm_8n14(x, y, z)) / ((16 - abs(-8))*(16 + abs(-8))); }
double Rm7n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm7n15(x, y, z) + (x*x + y*y + z*z)*Rm7n14(x, y, z)) / ((16 - abs(7))*(16 + abs(7))); }
double Rm_7n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_7n15(x, y, z) + (x*x + y*y + z*z)*Rm_7n14(x, y, z)) / ((16 - abs(-7))*(16 + abs(-7))); }
double Rm6n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm6n15(x, y, z) + (x*x + y*y + z*z)*Rm6n14(x, y, z)) / ((16 - abs(6))*(16 + abs(6))); }
double Rm_6n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_6n15(x, y, z) + (x*x + y*y + z*z)*Rm_6n14(x, y, z)) / ((16 - abs(-6))*(16 + abs(-6))); }
double Rm5n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm5n15(x, y, z) + (x*x + y*y + z*z)*Rm5n14(x, y, z)) / ((16 - abs(5))*(16 + abs(5))); }
double Rm_5n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_5n15(x, y, z) + (x*x + y*y + z*z)*Rm_5n14(x, y, z)) / ((16 - abs(-5))*(16 + abs(-5))); }
double Rm4n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm4n15(x, y, z) + (x*x + y*y + z*z)*Rm4n14(x, y, z)) / ((16 - abs(4))*(16 + abs(4))); }
double Rm_4n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_4n15(x, y, z) + (x*x + y*y + z*z)*Rm_4n14(x, y, z)) / ((16 - abs(-4))*(16 + abs(-4))); }
double Rm3n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm3n15(x, y, z) + (x*x + y*y + z*z)*Rm3n14(x, y, z)) / ((16 - abs(3))*(16 + abs(3))); }
double Rm_3n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_3n15(x, y, z) + (x*x + y*y + z*z)*Rm_3n14(x, y, z)) / ((16 - abs(-3))*(16 + abs(-3))); }
double Rm2n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm2n15(x, y, z) + (x*x + y*y + z*z)*Rm2n14(x, y, z)) / ((16 - abs(2))*(16 + abs(2))); }
double Rm_2n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_2n15(x, y, z) + (x*x + y*y + z*z)*Rm_2n14(x, y, z)) / ((16 - abs(-2))*(16 + abs(-2))); }
double Rm1n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm1n15(x, y, z) + (x*x + y*y + z*z)*Rm1n14(x, y, z)) / ((16 - abs(1))*(16 + abs(1))); }
double Rm_1n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm_1n15(x, y, z) + (x*x + y*y + z*z)*Rm_1n14(x, y, z)) / ((16 - abs(-1))*(16 + abs(-1))); }
double Rm0n16(double x, double y, double z) { return -((2 * 16 - 1)*z*Rm0n15(x, y, z) + (x*x + y*y + z*z)*Rm0n14(x, y, z)) / ((16 - abs(0))*(16 + abs(0))); }

double Rm17n17(double x, double y, double z) { return -(x*Rm16n16(x, y, z) + y*Rm_16n16(x, y, z)) / 34; }
double Rm_17n17(double x, double y, double z) { return (y*Rm16n16(x, y, z) - x*Rm_16n16(x, y, z)) / 34; }
double Rm16n17(double x, double y, double z) { return -z*Rm16n16(x, y, z); }
double Rm_16n17(double x, double y, double z) { return -z*Rm_16n16(x, y, z); }
double Rm15n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm15n16(x, y, z) + (x*x + y*y + z*z)*Rm15n15(x, y, z)) / ((17 - abs(15))*(17 + abs(15))); }
double Rm_15n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_15n16(x, y, z) + (x*x + y*y + z*z)*Rm_15n15(x, y, z)) / ((17 - abs(-15))*(17 + abs(-15))); }
double Rm14n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm14n16(x, y, z) + (x*x + y*y + z*z)*Rm14n15(x, y, z)) / ((17 - abs(14))*(17 + abs(14))); }
double Rm_14n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_14n16(x, y, z) + (x*x + y*y + z*z)*Rm_14n15(x, y, z)) / ((17 - abs(-14))*(17 + abs(-14))); }
double Rm13n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm13n16(x, y, z) + (x*x + y*y + z*z)*Rm13n15(x, y, z)) / ((17 - abs(13))*(17 + abs(13))); }
double Rm_13n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_13n16(x, y, z) + (x*x + y*y + z*z)*Rm_13n15(x, y, z)) / ((17 - abs(-13))*(17 + abs(-13))); }
double Rm12n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm12n16(x, y, z) + (x*x + y*y + z*z)*Rm12n15(x, y, z)) / ((17 - abs(12))*(17 + abs(12))); }
double Rm_12n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_12n16(x, y, z) + (x*x + y*y + z*z)*Rm_12n15(x, y, z)) / ((17 - abs(-12))*(17 + abs(-12))); }
double Rm11n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm11n16(x, y, z) + (x*x + y*y + z*z)*Rm11n15(x, y, z)) / ((17 - abs(11))*(17 + abs(11))); }
double Rm_11n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_11n16(x, y, z) + (x*x + y*y + z*z)*Rm_11n15(x, y, z)) / ((17 - abs(-11))*(17 + abs(-11))); }
double Rm10n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm10n16(x, y, z) + (x*x + y*y + z*z)*Rm10n15(x, y, z)) / ((17 - abs(10))*(17 + abs(10))); }
double Rm_10n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_10n16(x, y, z) + (x*x + y*y + z*z)*Rm_10n15(x, y, z)) / ((17 - abs(-10))*(17 + abs(-10))); }
double Rm9n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm9n16(x, y, z) + (x*x + y*y + z*z)*Rm9n15(x, y, z)) / ((17 - abs(9))*(17 + abs(9))); }
double Rm_9n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_9n16(x, y, z) + (x*x + y*y + z*z)*Rm_9n15(x, y, z)) / ((17 - abs(-9))*(17 + abs(-9))); }
double Rm8n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm8n16(x, y, z) + (x*x + y*y + z*z)*Rm8n15(x, y, z)) / ((17 - abs(8))*(17 + abs(8))); }
double Rm_8n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_8n16(x, y, z) + (x*x + y*y + z*z)*Rm_8n15(x, y, z)) / ((17 - abs(-8))*(17 + abs(-8))); }
double Rm7n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm7n16(x, y, z) + (x*x + y*y + z*z)*Rm7n15(x, y, z)) / ((17 - abs(7))*(17 + abs(7))); }
double Rm_7n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_7n16(x, y, z) + (x*x + y*y + z*z)*Rm_7n15(x, y, z)) / ((17 - abs(-7))*(17 + abs(-7))); }
double Rm6n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm6n16(x, y, z) + (x*x + y*y + z*z)*Rm6n15(x, y, z)) / ((17 - abs(6))*(17 + abs(6))); }
double Rm_6n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_6n16(x, y, z) + (x*x + y*y + z*z)*Rm_6n15(x, y, z)) / ((17 - abs(-6))*(17 + abs(-6))); }
double Rm5n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm5n16(x, y, z) + (x*x + y*y + z*z)*Rm5n15(x, y, z)) / ((17 - abs(5))*(17 + abs(5))); }
double Rm_5n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_5n16(x, y, z) + (x*x + y*y + z*z)*Rm_5n15(x, y, z)) / ((17 - abs(-5))*(17 + abs(-5))); }
double Rm4n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm4n16(x, y, z) + (x*x + y*y + z*z)*Rm4n15(x, y, z)) / ((17 - abs(4))*(17 + abs(4))); }
double Rm_4n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_4n16(x, y, z) + (x*x + y*y + z*z)*Rm_4n15(x, y, z)) / ((17 - abs(-4))*(17 + abs(-4))); }
double Rm3n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm3n16(x, y, z) + (x*x + y*y + z*z)*Rm3n15(x, y, z)) / ((17 - abs(3))*(17 + abs(3))); }
double Rm_3n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_3n16(x, y, z) + (x*x + y*y + z*z)*Rm_3n15(x, y, z)) / ((17 - abs(-3))*(17 + abs(-3))); }
double Rm2n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm2n16(x, y, z) + (x*x + y*y + z*z)*Rm2n15(x, y, z)) / ((17 - abs(2))*(17 + abs(2))); }
double Rm_2n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_2n16(x, y, z) + (x*x + y*y + z*z)*Rm_2n15(x, y, z)) / ((17 - abs(-2))*(17 + abs(-2))); }
double Rm1n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm1n16(x, y, z) + (x*x + y*y + z*z)*Rm1n15(x, y, z)) / ((17 - abs(1))*(17 + abs(1))); }
double Rm_1n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm_1n16(x, y, z) + (x*x + y*y + z*z)*Rm_1n15(x, y, z)) / ((17 - abs(-1))*(17 + abs(-1))); }
double Rm0n17(double x, double y, double z) { return -((2 * 17 - 1)*z*Rm0n16(x, y, z) + (x*x + y*y + z*z)*Rm0n15(x, y, z)) / ((17 - abs(0))*(17 + abs(0))); }

double Rm18n18(double x, double y, double z) { return -(x*Rm17n17(x, y, z) + y*Rm_17n17(x, y, z)) / 36; }
double Rm_18n18(double x, double y, double z) { return (y*Rm17n17(x, y, z) - x*Rm_17n17(x, y, z)) / 36; }
double Rm17n18(double x, double y, double z) { return -z*Rm17n17(x, y, z); }
double Rm_17n18(double x, double y, double z) { return -z*Rm_17n17(x, y, z); }
double Rm16n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm16n17(x, y, z) + (x*x + y*y + z*z)*Rm16n16(x, y, z)) / ((18 - abs(16))*(18 + abs(16))); }
double Rm_16n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_16n17(x, y, z) + (x*x + y*y + z*z)*Rm_16n16(x, y, z)) / ((18 - abs(-16))*(18 + abs(-16))); }
double Rm15n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm15n17(x, y, z) + (x*x + y*y + z*z)*Rm15n16(x, y, z)) / ((18 - abs(15))*(18 + abs(15))); }
double Rm_15n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_15n17(x, y, z) + (x*x + y*y + z*z)*Rm_15n16(x, y, z)) / ((18 - abs(-15))*(18 + abs(-15))); }
double Rm14n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm14n17(x, y, z) + (x*x + y*y + z*z)*Rm14n16(x, y, z)) / ((18 - abs(14))*(18 + abs(14))); }
double Rm_14n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_14n17(x, y, z) + (x*x + y*y + z*z)*Rm_14n16(x, y, z)) / ((18 - abs(-14))*(18 + abs(-14))); }
double Rm13n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm13n17(x, y, z) + (x*x + y*y + z*z)*Rm13n16(x, y, z)) / ((18 - abs(13))*(18 + abs(13))); }
double Rm_13n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_13n17(x, y, z) + (x*x + y*y + z*z)*Rm_13n16(x, y, z)) / ((18 - abs(-13))*(18 + abs(-13))); }
double Rm12n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm12n17(x, y, z) + (x*x + y*y + z*z)*Rm12n16(x, y, z)) / ((18 - abs(12))*(18 + abs(12))); }
double Rm_12n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_12n17(x, y, z) + (x*x + y*y + z*z)*Rm_12n16(x, y, z)) / ((18 - abs(-12))*(18 + abs(-12))); }
double Rm11n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm11n17(x, y, z) + (x*x + y*y + z*z)*Rm11n16(x, y, z)) / ((18 - abs(11))*(18 + abs(11))); }
double Rm_11n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_11n17(x, y, z) + (x*x + y*y + z*z)*Rm_11n16(x, y, z)) / ((18 - abs(-11))*(18 + abs(-11))); }
double Rm10n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm10n17(x, y, z) + (x*x + y*y + z*z)*Rm10n16(x, y, z)) / ((18 - abs(10))*(18 + abs(10))); }
double Rm_10n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_10n17(x, y, z) + (x*x + y*y + z*z)*Rm_10n16(x, y, z)) / ((18 - abs(-10))*(18 + abs(-10))); }
double Rm9n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm9n17(x, y, z) + (x*x + y*y + z*z)*Rm9n16(x, y, z)) / ((18 - abs(9))*(18 + abs(9))); }
double Rm_9n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_9n17(x, y, z) + (x*x + y*y + z*z)*Rm_9n16(x, y, z)) / ((18 - abs(-9))*(18 + abs(-9))); }
double Rm8n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm8n17(x, y, z) + (x*x + y*y + z*z)*Rm8n16(x, y, z)) / ((18 - abs(8))*(18 + abs(8))); }
double Rm_8n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_8n17(x, y, z) + (x*x + y*y + z*z)*Rm_8n16(x, y, z)) / ((18 - abs(-8))*(18 + abs(-8))); }
double Rm7n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm7n17(x, y, z) + (x*x + y*y + z*z)*Rm7n16(x, y, z)) / ((18 - abs(7))*(18 + abs(7))); }
double Rm_7n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_7n17(x, y, z) + (x*x + y*y + z*z)*Rm_7n16(x, y, z)) / ((18 - abs(-7))*(18 + abs(-7))); }
double Rm6n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm6n17(x, y, z) + (x*x + y*y + z*z)*Rm6n16(x, y, z)) / ((18 - abs(6))*(18 + abs(6))); }
double Rm_6n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_6n17(x, y, z) + (x*x + y*y + z*z)*Rm_6n16(x, y, z)) / ((18 - abs(-6))*(18 + abs(-6))); }
double Rm5n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm5n17(x, y, z) + (x*x + y*y + z*z)*Rm5n16(x, y, z)) / ((18 - abs(5))*(18 + abs(5))); }
double Rm_5n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_5n17(x, y, z) + (x*x + y*y + z*z)*Rm_5n16(x, y, z)) / ((18 - abs(-5))*(18 + abs(-5))); }
double Rm4n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm4n17(x, y, z) + (x*x + y*y + z*z)*Rm4n16(x, y, z)) / ((18 - abs(4))*(18 + abs(4))); }
double Rm_4n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_4n17(x, y, z) + (x*x + y*y + z*z)*Rm_4n16(x, y, z)) / ((18 - abs(-4))*(18 + abs(-4))); }
double Rm3n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm3n17(x, y, z) + (x*x + y*y + z*z)*Rm3n16(x, y, z)) / ((18 - abs(3))*(18 + abs(3))); }
double Rm_3n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_3n17(x, y, z) + (x*x + y*y + z*z)*Rm_3n16(x, y, z)) / ((18 - abs(-3))*(18 + abs(-3))); }
double Rm2n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm2n17(x, y, z) + (x*x + y*y + z*z)*Rm2n16(x, y, z)) / ((18 - abs(2))*(18 + abs(2))); }
double Rm_2n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_2n17(x, y, z) + (x*x + y*y + z*z)*Rm_2n16(x, y, z)) / ((18 - abs(-2))*(18 + abs(-2))); }
double Rm1n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm1n17(x, y, z) + (x*x + y*y + z*z)*Rm1n16(x, y, z)) / ((18 - abs(1))*(18 + abs(1))); }
double Rm_1n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm_1n17(x, y, z) + (x*x + y*y + z*z)*Rm_1n16(x, y, z)) / ((18 - abs(-1))*(18 + abs(-1))); }
double Rm0n18(double x, double y, double z) { return -((2 * 18 - 1)*z*Rm0n17(x, y, z) + (x*x + y*y + z*z)*Rm0n16(x, y, z)) / ((18 - abs(0))*(18 + abs(0))); }

double Rm19n19(double x, double y, double z) { return -(x*Rm18n18(x, y, z) + y*Rm_18n18(x, y, z)) / 38; }
double Rm_19n19(double x, double y, double z) { return (y*Rm18n18(x, y, z) - x*Rm_18n18(x, y, z)) / 38; }
double Rm18n19(double x, double y, double z) { return -z*Rm18n18(x, y, z); }
double Rm_18n19(double x, double y, double z) { return -z*Rm_18n18(x, y, z); }
double Rm17n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm17n18(x, y, z) + (x*x + y*y + z*z)*Rm17n17(x, y, z)) / ((19 - abs(17))*(19 + abs(17))); }
double Rm_17n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_17n18(x, y, z) + (x*x + y*y + z*z)*Rm_17n17(x, y, z)) / ((19 - abs(-17))*(19 + abs(-17))); }
double Rm16n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm16n18(x, y, z) + (x*x + y*y + z*z)*Rm16n17(x, y, z)) / ((19 - abs(16))*(19 + abs(16))); }
double Rm_16n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_16n18(x, y, z) + (x*x + y*y + z*z)*Rm_16n17(x, y, z)) / ((19 - abs(-16))*(19 + abs(-16))); }
double Rm15n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm15n18(x, y, z) + (x*x + y*y + z*z)*Rm15n17(x, y, z)) / ((19 - abs(15))*(19 + abs(15))); }
double Rm_15n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_15n18(x, y, z) + (x*x + y*y + z*z)*Rm_15n17(x, y, z)) / ((19 - abs(-15))*(19 + abs(-15))); }
double Rm14n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm14n18(x, y, z) + (x*x + y*y + z*z)*Rm14n17(x, y, z)) / ((19 - abs(14))*(19 + abs(14))); }
double Rm_14n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_14n18(x, y, z) + (x*x + y*y + z*z)*Rm_14n17(x, y, z)) / ((19 - abs(-14))*(19 + abs(-14))); }
double Rm13n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm13n18(x, y, z) + (x*x + y*y + z*z)*Rm13n17(x, y, z)) / ((19 - abs(13))*(19 + abs(13))); }
double Rm_13n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_13n18(x, y, z) + (x*x + y*y + z*z)*Rm_13n17(x, y, z)) / ((19 - abs(-13))*(19 + abs(-13))); }
double Rm12n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm12n18(x, y, z) + (x*x + y*y + z*z)*Rm12n17(x, y, z)) / ((19 - abs(12))*(19 + abs(12))); }
double Rm_12n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_12n18(x, y, z) + (x*x + y*y + z*z)*Rm_12n17(x, y, z)) / ((19 - abs(-12))*(19 + abs(-12))); }
double Rm11n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm11n18(x, y, z) + (x*x + y*y + z*z)*Rm11n17(x, y, z)) / ((19 - abs(11))*(19 + abs(11))); }
double Rm_11n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_11n18(x, y, z) + (x*x + y*y + z*z)*Rm_11n17(x, y, z)) / ((19 - abs(-11))*(19 + abs(-11))); }
double Rm10n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm10n18(x, y, z) + (x*x + y*y + z*z)*Rm10n17(x, y, z)) / ((19 - abs(10))*(19 + abs(10))); }
double Rm_10n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_10n18(x, y, z) + (x*x + y*y + z*z)*Rm_10n17(x, y, z)) / ((19 - abs(-10))*(19 + abs(-10))); }
double Rm9n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm9n18(x, y, z) + (x*x + y*y + z*z)*Rm9n17(x, y, z)) / ((19 - abs(9))*(19 + abs(9))); }
double Rm_9n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_9n18(x, y, z) + (x*x + y*y + z*z)*Rm_9n17(x, y, z)) / ((19 - abs(-9))*(19 + abs(-9))); }
double Rm8n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm8n18(x, y, z) + (x*x + y*y + z*z)*Rm8n17(x, y, z)) / ((19 - abs(8))*(19 + abs(8))); }
double Rm_8n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_8n18(x, y, z) + (x*x + y*y + z*z)*Rm_8n17(x, y, z)) / ((19 - abs(-8))*(19 + abs(-8))); }
double Rm7n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm7n18(x, y, z) + (x*x + y*y + z*z)*Rm7n17(x, y, z)) / ((19 - abs(7))*(19 + abs(7))); }
double Rm_7n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_7n18(x, y, z) + (x*x + y*y + z*z)*Rm_7n17(x, y, z)) / ((19 - abs(-7))*(19 + abs(-7))); }
double Rm6n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm6n18(x, y, z) + (x*x + y*y + z*z)*Rm6n17(x, y, z)) / ((19 - abs(6))*(19 + abs(6))); }
double Rm_6n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_6n18(x, y, z) + (x*x + y*y + z*z)*Rm_6n17(x, y, z)) / ((19 - abs(-6))*(19 + abs(-6))); }
double Rm5n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm5n18(x, y, z) + (x*x + y*y + z*z)*Rm5n17(x, y, z)) / ((19 - abs(5))*(19 + abs(5))); }
double Rm_5n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_5n18(x, y, z) + (x*x + y*y + z*z)*Rm_5n17(x, y, z)) / ((19 - abs(-5))*(19 + abs(-5))); }
double Rm4n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm4n18(x, y, z) + (x*x + y*y + z*z)*Rm4n17(x, y, z)) / ((19 - abs(4))*(19 + abs(4))); }
double Rm_4n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_4n18(x, y, z) + (x*x + y*y + z*z)*Rm_4n17(x, y, z)) / ((19 - abs(-4))*(19 + abs(-4))); }
double Rm3n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm3n18(x, y, z) + (x*x + y*y + z*z)*Rm3n17(x, y, z)) / ((19 - abs(3))*(19 + abs(3))); }
double Rm_3n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_3n18(x, y, z) + (x*x + y*y + z*z)*Rm_3n17(x, y, z)) / ((19 - abs(-3))*(19 + abs(-3))); }
double Rm2n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm2n18(x, y, z) + (x*x + y*y + z*z)*Rm2n17(x, y, z)) / ((19 - abs(2))*(19 + abs(2))); }
double Rm_2n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_2n18(x, y, z) + (x*x + y*y + z*z)*Rm_2n17(x, y, z)) / ((19 - abs(-2))*(19 + abs(-2))); }
double Rm1n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm1n18(x, y, z) + (x*x + y*y + z*z)*Rm1n17(x, y, z)) / ((19 - abs(1))*(19 + abs(1))); }
double Rm_1n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm_1n18(x, y, z) + (x*x + y*y + z*z)*Rm_1n17(x, y, z)) / ((19 - abs(-1))*(19 + abs(-1))); }
double Rm0n19(double x, double y, double z) { return -((2 * 19 - 1)*z*Rm0n18(x, y, z) + (x*x + y*y + z*z)*Rm0n17(x, y, z)) / ((19 - abs(0))*(19 + abs(0))); }

//set R^{heat}mn  Normalization
double RRm0n0(double x, double y, double z) { return sqrt(fac(0)*fac(0))*Rm0n0(x,y,z); }

double RRm1n1(double x, double y, double z) { return sqrt(fac(0)*fac(2))*Rm1n1(x, y, z); }
double RRm_1n1(double x, double y, double z) { return sqrt(fac(2)*fac(0))*Rm_1n1(x, y, z); }
double RRm0n1(double x, double y, double z) { return sqrt(fac(1)*fac(1))*Rm0n1(x, y, z); }

double RRm2n2(double x, double y, double z) { return sqrt(fac(0)*fac(4))*Rm2n2(x, y, z); }
double RRm_2n2(double x, double y, double z) { return sqrt(fac(4)*fac(0))*Rm_2n2(x, y, z); }
double RRm1n2(double x, double y, double z) { return sqrt(fac(1)*fac(3))*Rm1n2(x, y, z); }
double RRm_1n2(double x, double y, double z) { return sqrt(fac(3)*fac(1))*Rm_1n2(x, y, z); }
double RRm0n2(double x, double y, double z) { return sqrt(fac(2)*fac(2))*Rm0n2(x, y, z); }

double RRm3n3(double x, double y, double z) { return sqrt(fac(0)*fac(6))*Rm3n3(x, y, z); }
double RRm_3n3(double x, double y, double z) { return sqrt(fac(6)*fac(0))*Rm_3n3(x, y, z); }
double RRm2n3(double x, double y, double z) { return sqrt(fac(1)*fac(5))*Rm2n3(x, y, z); }
double RRm_2n3(double x, double y, double z) { return sqrt(fac(5)*fac(1))*Rm_2n3(x, y, z); }
double RRm1n3(double x, double y, double z) { return sqrt(fac(2)*fac(4))*Rm1n3(x, y, z); }
double RRm_1n3(double x, double y, double z) { return sqrt(fac(4)*fac(2))*Rm_1n3(x, y, z); }
double RRm0n3(double x, double y, double z) { return sqrt(fac(3)*fac(3))*Rm0n3(x, y, z); }

double RRm4n4(double x, double y, double z) { return sqrt(fac(0)*fac(8))*Rm4n4(x, y, z); }
double RRm_4n4(double x, double y, double z) { return sqrt(fac(8)*fac(0))*Rm_4n4(x, y, z); }
double RRm3n4(double x, double y, double z) { return sqrt(fac(1)*fac(7))*Rm3n4(x, y, z); }
double RRm_3n4(double x, double y, double z) { return sqrt(fac(7)*fac(1))*Rm_3n4(x, y, z); }
double RRm2n4(double x, double y, double z) { return sqrt(fac(2)*fac(6))*Rm2n4(x, y, z); }
double RRm_2n4(double x, double y, double z) { return sqrt(fac(6)*fac(2))*Rm_2n4(x, y, z); }
double RRm1n4(double x, double y, double z) { return sqrt(fac(3)*fac(5))*Rm1n4(x, y, z); }
double RRm_1n4(double x, double y, double z) { return sqrt(fac(5)*fac(3))*Rm_1n4(x, y, z); }
double RRm0n4(double x, double y, double z) { return sqrt(fac(4)*fac(4))*Rm0n4(x, y, z); }

double RRm5n5(double x, double y, double z) { return sqrt(fac(0)*fac(10))*Rm5n5(x, y, z); }
double RRm_5n5(double x, double y, double z) { return sqrt(fac(10)*fac(0))*Rm_5n5(x, y, z); }
double RRm4n5(double x, double y, double z) { return sqrt(fac(1)*fac(9))*Rm4n5(x, y, z); }
double RRm_4n5(double x, double y, double z) { return sqrt(fac(9)*fac(1))*Rm_4n5(x, y, z); }
double RRm3n5(double x, double y, double z) { return sqrt(fac(2)*fac(8))*Rm3n5(x, y, z); }
double RRm_3n5(double x, double y, double z) { return sqrt(fac(8)*fac(2))*Rm_3n5(x, y, z); }
double RRm2n5(double x, double y, double z) { return sqrt(fac(3)*fac(7))*Rm2n5(x, y, z); }
double RRm_2n5(double x, double y, double z) { return sqrt(fac(7)*fac(3))*Rm_2n5(x, y, z); }
double RRm1n5(double x, double y, double z) { return sqrt(fac(4)*fac(6))*Rm1n5(x, y, z); }
double RRm_1n5(double x, double y, double z) { return sqrt(fac(6)*fac(4))*Rm_1n5(x, y, z); }
double RRm0n5(double x, double y, double z) { return sqrt(fac(5)*fac(5))*Rm0n5(x, y, z); }

double RRm6n6(double x, double y, double z) { return sqrt(fac(0)*fac(12))*Rm6n6(x, y, z); }
double RRm_6n6(double x, double y, double z) { return sqrt(fac(12)*fac(0))*Rm_6n6(x, y, z); }
double RRm5n6(double x, double y, double z) { return sqrt(fac(1)*fac(11))*Rm5n6(x, y, z); }
double RRm_5n6(double x, double y, double z) { return sqrt(fac(11)*fac(1))*Rm_5n6(x, y, z); }
double RRm4n6(double x, double y, double z) { return sqrt(fac(2)*fac(10))*Rm4n6(x, y, z); }
double RRm_4n6(double x, double y, double z) { return sqrt(fac(10)*fac(2))*Rm_4n6(x, y, z); }
double RRm3n6(double x, double y, double z) { return sqrt(fac(3)*fac(9))*Rm3n6(x, y, z); }
double RRm_3n6(double x, double y, double z) { return sqrt(fac(9)*fac(3))*Rm_3n6(x, y, z); }
double RRm2n6(double x, double y, double z) { return sqrt(fac(4)*fac(8))*Rm2n6(x, y, z); }
double RRm_2n6(double x, double y, double z) { return sqrt(fac(8)*fac(4))*Rm_2n6(x, y, z); }
double RRm1n6(double x, double y, double z) { return sqrt(fac(5)*fac(7))*Rm1n6(x, y, z); }
double RRm_1n6(double x, double y, double z) { return sqrt(fac(7)*fac(5))*Rm_1n6(x, y, z); }
double RRm0n6(double x, double y, double z) { return sqrt(fac(6)*fac(6))*Rm0n6(x, y, z); }

double RRm7n7(double x, double y, double z) { return sqrt(fac(0)*fac(14))*Rm7n7(x, y, z); }
double RRm_7n7(double x, double y, double z) { return sqrt(fac(14)*fac(0))*Rm_7n7(x, y, z); }
double RRm6n7(double x, double y, double z) { return sqrt(fac(1)*fac(13))*Rm6n7(x, y, z); }
double RRm_6n7(double x, double y, double z) { return sqrt(fac(13)*fac(1))*Rm_6n7(x, y, z); }
double RRm5n7(double x, double y, double z) { return sqrt(fac(2)*fac(12))*Rm5n7(x, y, z); }
double RRm_5n7(double x, double y, double z) { return sqrt(fac(12)*fac(2))*Rm_5n7(x, y, z); }
double RRm4n7(double x, double y, double z) { return sqrt(fac(3)*fac(11))*Rm4n7(x, y, z); }
double RRm_4n7(double x, double y, double z) { return sqrt(fac(11)*fac(3))*Rm_4n7(x, y, z); }
double RRm3n7(double x, double y, double z) { return sqrt(fac(4)*fac(10))*Rm3n7(x, y, z); }
double RRm_3n7(double x, double y, double z) { return sqrt(fac(10)*fac(4))*Rm_3n7(x, y, z); }
double RRm2n7(double x, double y, double z) { return sqrt(fac(5)*fac(9))*Rm2n7(x, y, z); }
double RRm_2n7(double x, double y, double z) { return sqrt(fac(9)*fac(5))*Rm_2n7(x, y, z); }
double RRm1n7(double x, double y, double z) { return sqrt(fac(6)*fac(8))*Rm1n7(x, y, z); }
double RRm_1n7(double x, double y, double z) { return sqrt(fac(8)*fac(6))*Rm_1n7(x, y, z); }
double RRm0n7(double x, double y, double z) { return sqrt(fac(7)*fac(7))*Rm0n7(x, y, z); }

double RRm8n8(double x, double y, double z) { return sqrt(fac(0.0)*fac(16.0))*Rm8n8(x, y, z); }
double RRm_8n8(double x, double y, double z) { return sqrt(fac(16.0)*fac(0.0))*Rm_8n8(x, y, z); }
double RRm7n8(double x, double y, double z) { return sqrt(fac(1.0)*fac(15.0))*Rm7n8(x, y, z); }
double RRm_7n8(double x, double y, double z) { return sqrt(fac(15.0)*fac(1.0))*Rm_7n8(x, y, z); }
double RRm6n8(double x, double y, double z) {return sqrt(fac(2.0)*fac(14.0))*Rm6n8(x, y, z); }
double RRm_6n8(double x, double y, double z) { return sqrt(fac(14.0)*fac(2.0))*Rm_6n8(x, y, z); }
double RRm5n8(double x, double y, double z) { return sqrt(fac(3.0)*fac(13.0))*Rm5n8(x, y, z); }
double RRm_5n8(double x, double y, double z) { return sqrt(fac(13.0)*fac(3.0))*Rm_5n8(x, y, z); }
double RRm4n8(double x, double y, double z) { return sqrt(fac(4.0)*fac(12.0))*Rm4n8(x, y, z); }
double RRm_4n8(double x, double y, double z) { return sqrt(fac(12.0)*fac(4.0))*Rm_4n8(x, y, z); }
double RRm3n8(double x, double y, double z) { return sqrt(fac(5.0)*fac(11.0))*Rm3n8(x, y, z); }
double RRm_3n8(double x, double y, double z) { return sqrt(fac(11.0)*fac(5.0))*Rm_3n8(x, y, z); }
double RRm2n8(double x, double y, double z) { return sqrt(fac(6.0)*fac(10.0))*Rm2n8(x, y, z); }
double RRm_2n8(double x, double y, double z) { return sqrt(fac(10.0)*fac(6.0))*Rm_2n8(x, y, z); }
double RRm1n8(double x, double y, double z) { return sqrt(fac(7.0)*fac(9.0))*Rm1n8(x, y, z); }
double RRm_1n8(double x, double y, double z) { return sqrt(fac(9.0)*fac(7.0))*Rm_1n8(x, y, z); }
double RRm0n8(double x, double y, double z) { return sqrt(fac(8.0)*fac(8.0))*Rm0n8(x, y, z); }

double RRm9n9(double x, double y, double z) { return sqrt(fac(0.0)*fac(18.0))*Rm9n9(x, y, z); }
double RRm_9n9(double x, double y, double z) { return sqrt(fac(18.0)*fac(0.0))*Rm_9n9(x, y, z); }
double RRm8n9(double x, double y, double z) { return sqrt(fac(1.0)*fac(17.0))*Rm8n9(x, y, z); }
double RRm_8n9(double x, double y, double z) { return sqrt(fac(17.0)*fac(1.0))*Rm_8n9(x, y, z); }
double RRm7n9(double x, double y, double z) { return sqrt(fac(2.0)*fac(16.0))*Rm7n9(x, y, z); }
double RRm_7n9(double x, double y, double z) { return sqrt(fac(16.0)*fac(2.0))*Rm_7n9(x, y, z); }
double RRm6n9(double x, double y, double z) { return sqrt(fac(3.0)*fac(15.0))*Rm6n9(x, y, z); }
double RRm_6n9(double x, double y, double z) { return sqrt(fac(15.0)*fac(3.0))*Rm_6n9(x, y, z); }
double RRm5n9(double x, double y, double z) { return sqrt(fac(4.0)*fac(14.0))*Rm5n9(x, y, z); }
double RRm_5n9(double x, double y, double z) { return sqrt(fac(14.0)*fac(4.0))*Rm_5n9(x, y, z); }
double RRm4n9(double x, double y, double z) { return sqrt(fac(5.0)*fac(13.0))*Rm4n9(x, y, z); }
double RRm_4n9(double x, double y, double z) { return sqrt(fac(13.0)*fac(5.0))*Rm_4n9(x, y, z); }
double RRm3n9(double x, double y, double z) { return sqrt(fac(6.0)*fac(12.0))*Rm3n9(x, y, z); }
double RRm_3n9(double x, double y, double z) { return sqrt(fac(12.0)*fac(6.0))*Rm_3n9(x, y, z); }
double RRm2n9(double x, double y, double z) { return sqrt(fac(7.0)*fac(11.0))*Rm2n9(x, y, z); }
double RRm_2n9(double x, double y, double z) { return sqrt(fac(11.0)*fac(7.0))*Rm_2n9(x, y, z); }
double RRm1n9(double x, double y, double z) { return sqrt(fac(8.0)*fac(10.0))*Rm1n9(x, y, z); }
double RRm_1n9(double x, double y, double z) { return sqrt(fac(10.0)*fac(8.0))*Rm_1n9(x, y, z); }
double RRm0n9(double x, double y, double z) { return sqrt(fac(9.0)*fac(9.0))*Rm0n9(x, y, z); }

double RRm10n10(double x, double y, double z) { return sqrt(fac(0.0)*fac(20.0))*Rm10n10(x, y, z); }
double RRm_10n10(double x, double y, double z) { return sqrt(fac(20.0)*fac(0.0))*Rm_10n10(x, y, z); }
double RRm9n10(double x, double y, double z) { return sqrt(fac(1.0)*fac(19.0))*Rm9n10(x, y, z); }
double RRm_9n10(double x, double y, double z) { return sqrt(fac(19.0)*fac(1.0))*Rm_9n10(x, y, z); }
double RRm8n10(double x, double y, double z) { return sqrt(fac(2.0)*fac(18.0))*Rm8n10(x, y, z); }
double RRm_8n10(double x, double y, double z) { return sqrt(fac(18.0)*fac(2.0))*Rm_8n10(x, y, z); }
double RRm7n10(double x, double y, double z) { return sqrt(fac(3.0)*fac(17.0))*Rm7n10(x, y, z); }
double RRm_7n10(double x, double y, double z) { return sqrt(fac(17.0)*fac(3.0))*Rm_7n10(x, y, z); }
double RRm6n10(double x, double y, double z) { return sqrt(fac(4.0)*fac(16.0))*Rm6n10(x, y, z); }
double RRm_6n10(double x, double y, double z) { return sqrt(fac(16.0)*fac(4.0))*Rm_6n10(x, y, z); }
double RRm5n10(double x, double y, double z) { return sqrt(fac(5.0)*fac(15.0))*Rm5n10(x, y, z); }
double RRm_5n10(double x, double y, double z) { return sqrt(fac(15.0)*fac(5.0))*Rm_5n10(x, y, z); }
double RRm4n10(double x, double y, double z) { return sqrt(fac(6.0)*fac(14.0))*Rm4n10(x, y, z); }
double RRm_4n10(double x, double y, double z) { return sqrt(fac(14.0)*fac(6.0))*Rm_4n10(x, y, z); }
double RRm3n10(double x, double y, double z) { return sqrt(fac(7.0)*fac(13.0))*Rm3n10(x, y, z); }
double RRm_3n10(double x, double y, double z) { return sqrt(fac(13.0)*fac(7.0))*Rm_3n10(x, y, z); }
double RRm2n10(double x, double y, double z) { return sqrt(fac(8.0)*fac(12.0))*Rm2n10(x, y, z); }
double RRm_2n10(double x, double y, double z) { return sqrt(fac(12.0)*fac(8.0))*Rm_2n10(x, y, z); }
double RRm1n10(double x, double y, double z) { return sqrt(fac(9.0)*fac(11.0))*Rm1n10(x, y, z); }
double RRm_1n10(double x, double y, double z) { return sqrt(fac(11.0)*fac(9.0))*Rm_1n10(x, y, z); }
double RRm0n10(double x, double y, double z) { return sqrt(fac(10.0)*fac(10.0))*Rm0n10(x, y, z); }

double RRm11n11(double x, double y, double z) { return sqrt(fac(0.0)*fac(22.0))*Rm11n11(x, y, z); }
double RRm_11n11(double x, double y, double z) { return sqrt(fac(22.0)*fac(0.0))*Rm_11n11(x, y, z); }
double RRm10n11(double x, double y, double z) { return sqrt(fac(1.0)*fac(21.0))*Rm10n11(x, y, z); }
double RRm_10n11(double x, double y, double z) { return sqrt(fac(21.0)*fac(1.0))*Rm_10n11(x, y, z); }
double RRm9n11(double x, double y, double z) { return sqrt(fac(2.0)*fac(20.0))*Rm9n11(x, y, z); }
double RRm_9n11(double x, double y, double z) { return sqrt(fac(20.0)*fac(2.0))*Rm_9n11(x, y, z); }
double RRm8n11(double x, double y, double z) { return sqrt(fac(3.0)*fac(19.0))*Rm8n11(x, y, z); }
double RRm_8n11(double x, double y, double z) { return sqrt(fac(19.0)*fac(3.0))*Rm_8n11(x, y, z); }
double RRm7n11(double x, double y, double z) { return sqrt(fac(4.0)*fac(18.0))*Rm7n11(x, y, z); }
double RRm_7n11(double x, double y, double z) { return sqrt(fac(18.0)*fac(4.0))*Rm_7n11(x, y, z); }
double RRm6n11(double x, double y, double z) { return sqrt(fac(5.0)*fac(17.0))*Rm6n11(x, y, z); }
double RRm_6n11(double x, double y, double z) { return sqrt(fac(17.0)*fac(5.0))*Rm_6n11(x, y, z); }
double RRm5n11(double x, double y, double z) { return sqrt(fac(6.0)*fac(16.0))*Rm5n11(x, y, z); }
double RRm_5n11(double x, double y, double z) { return sqrt(fac(16.0)*fac(6.0))*Rm_5n11(x, y, z); }
double RRm4n11(double x, double y, double z) { return sqrt(fac(7.0)*fac(15.0))*Rm4n11(x, y, z); }
double RRm_4n11(double x, double y, double z) { return sqrt(fac(15.0)*fac(7.0))*Rm_4n11(x, y, z); }
double RRm3n11(double x, double y, double z) { return sqrt(fac(8.0)*fac(14.0))*Rm3n11(x, y, z); }
double RRm_3n11(double x, double y, double z) { return sqrt(fac(14.0)*fac(8.0))*Rm_3n11(x, y, z); }
double RRm2n11(double x, double y, double z) { return sqrt(fac(9.0)*fac(13.0))*Rm2n11(x, y, z); }
double RRm_2n11(double x, double y, double z) { return sqrt(fac(13.0)*fac(9.0))*Rm_2n11(x, y, z); }
double RRm1n11(double x, double y, double z) { return sqrt(fac(10.0)*fac(12.0))*Rm1n11(x, y, z); }
double RRm_1n11(double x, double y, double z) { return sqrt(fac(12.0)*fac(10.0))*Rm_1n11(x, y, z); }
double RRm0n11(double x, double y, double z) { return sqrt(fac(11.0)*fac(11.0))*Rm0n11(x, y, z); }

double RRm12n12(double x, double y, double z) { return sqrt(fac(0.0)*fac(24.0))*Rm12n12(x, y, z); }
double RRm_12n12(double x, double y, double z) { return sqrt(fac(24.0)*fac(0.0))*Rm_12n12(x, y, z); }
double RRm11n12(double x, double y, double z) { return sqrt(fac(1.0)*fac(23.0))*Rm11n12(x, y, z); }
double RRm_11n12(double x, double y, double z) { return sqrt(fac(23.0)*fac(1.0))*Rm_11n12(x, y, z); }
double RRm10n12(double x, double y, double z) { return sqrt(fac(2.0)*fac(22.0))*Rm10n12(x, y, z); }
double RRm_10n12(double x, double y, double z) { return sqrt(fac(22.0)*fac(2.0))*Rm_10n12(x, y, z); }
double RRm9n12(double x, double y, double z) { return sqrt(fac(3.0)*fac(21.0))*Rm9n12(x, y, z); }
double RRm_9n12(double x, double y, double z) { return sqrt(fac(21.0)*fac(3.0))*Rm_9n12(x, y, z); }
double RRm8n12(double x, double y, double z) { return sqrt(fac(4.0)*fac(20.0))*Rm8n12(x, y, z); }
double RRm_8n12(double x, double y, double z) { return sqrt(fac(20.0)*fac(4.0))*Rm_8n12(x, y, z); }
double RRm7n12(double x, double y, double z) { return sqrt(fac(5.0)*fac(19.0))*Rm7n12(x, y, z); }
double RRm_7n12(double x, double y, double z) { return sqrt(fac(19.0)*fac(5.0))*Rm_7n12(x, y, z); }
double RRm6n12(double x, double y, double z) { return sqrt(fac(6.0)*fac(18.0))*Rm6n12(x, y, z); }
double RRm_6n12(double x, double y, double z) { return sqrt(fac(18.0)*fac(6.0))*Rm_6n12(x, y, z); }
double RRm5n12(double x, double y, double z) { return sqrt(fac(7.0)*fac(17.0))*Rm5n12(x, y, z); }
double RRm_5n12(double x, double y, double z) { return sqrt(fac(17.0)*fac(7.0))*Rm_5n12(x, y, z); }
double RRm4n12(double x, double y, double z) { return sqrt(fac(8.0)*fac(16.0))*Rm4n12(x, y, z); }
double RRm_4n12(double x, double y, double z) { return sqrt(fac(16.0)*fac(8.0))*Rm_4n12(x, y, z); }
double RRm3n12(double x, double y, double z) { return sqrt(fac(9.0)*fac(15.0))*Rm3n12(x, y, z); }
double RRm_3n12(double x, double y, double z) { return sqrt(fac(15.0)*fac(9.0))*Rm_3n12(x, y, z); }
double RRm2n12(double x, double y, double z) { return sqrt(fac(10.0)*fac(14.0))*Rm2n12(x, y, z); }
double RRm_2n12(double x, double y, double z) { return sqrt(fac(14.0)*fac(10.0))*Rm_2n12(x, y, z); }
double RRm1n12(double x, double y, double z) { return sqrt(fac(11.0)*fac(13.0))*Rm1n12(x, y, z); }
double RRm_1n12(double x, double y, double z) { return sqrt(fac(13.0)*fac(11.0))*Rm_1n12(x, y, z); }
double RRm0n12(double x, double y, double z) { return sqrt(fac(12.0)*fac(12.0))*Rm0n12(x, y, z); }

double RRm13n13(double x, double y, double z) { return sqrt(fac(0.0)*fac(26.0))*Rm13n13(x, y, z); }
double RRm_13n13(double x, double y, double z) { return sqrt(fac(26.0)*fac(0.0))*Rm_13n13(x, y, z); }
double RRm12n13(double x, double y, double z) { return sqrt(fac(1.0)*fac(25.0))*Rm12n13(x, y, z); }
double RRm_12n13(double x, double y, double z) { return sqrt(fac(25.0)*fac(1.0))*Rm_12n13(x, y, z); }
double RRm11n13(double x, double y, double z) { return sqrt(fac(2.0)*fac(24.0))*Rm11n13(x, y, z); }
double RRm_11n13(double x, double y, double z) { return sqrt(fac(24.0)*fac(2.0))*Rm_11n13(x, y, z); }
double RRm10n13(double x, double y, double z) { return sqrt(fac(3.0)*fac(23.0))*Rm10n13(x, y, z); }
double RRm_10n13(double x, double y, double z) { return sqrt(fac(23.0)*fac(3.0))*Rm_10n13(x, y, z); }
double RRm9n13(double x, double y, double z) { return sqrt(fac(4.0)*fac(22.0))*Rm9n13(x, y, z); }
double RRm_9n13(double x, double y, double z) { return sqrt(fac(22.0)*fac(4.0))*Rm_9n13(x, y, z); }
double RRm8n13(double x, double y, double z) { return sqrt(fac(5.0)*fac(21.0))*Rm8n13(x, y, z); }
double RRm_8n13(double x, double y, double z) { return sqrt(fac(21.0)*fac(5.0))*Rm_8n13(x, y, z); }
double RRm7n13(double x, double y, double z) { return sqrt(fac(6.0)*fac(20.0))*Rm7n13(x, y, z); }
double RRm_7n13(double x, double y, double z) { return sqrt(fac(20.0)*fac(6.0))*Rm_7n13(x, y, z); }
double RRm6n13(double x, double y, double z) { return sqrt(fac(7.0)*fac(19.0))*Rm6n13(x, y, z); }
double RRm_6n13(double x, double y, double z) { return sqrt(fac(19.0)*fac(7.0))*Rm_6n13(x, y, z); }
double RRm5n13(double x, double y, double z) { return sqrt(fac(8.0)*fac(18.0))*Rm5n13(x, y, z); }
double RRm_5n13(double x, double y, double z) { return sqrt(fac(18.0)*fac(8.0))*Rm_5n13(x, y, z); }
double RRm4n13(double x, double y, double z) { return sqrt(fac(9.0)*fac(17.0))*Rm4n13(x, y, z); }
double RRm_4n13(double x, double y, double z) { return sqrt(fac(17.0)*fac(9.0))*Rm_4n13(x, y, z); }
double RRm3n13(double x, double y, double z) { return sqrt(fac(10.0)*fac(16.0))*Rm3n13(x, y, z); }
double RRm_3n13(double x, double y, double z) { return sqrt(fac(16.0)*fac(10.0))*Rm_3n13(x, y, z); }
double RRm2n13(double x, double y, double z) { return sqrt(fac(11.0)*fac(15.0))*Rm2n13(x, y, z); }
double RRm_2n13(double x, double y, double z) { return sqrt(fac(15.0)*fac(11.0))*Rm_2n13(x, y, z); }
double RRm1n13(double x, double y, double z) { return sqrt(fac(12.0)*fac(14.0))*Rm1n13(x, y, z); }
double RRm_1n13(double x, double y, double z) { return sqrt(fac(14.0)*fac(12.0))*Rm_1n13(x, y, z); }
double RRm0n13(double x, double y, double z) { return sqrt(fac(13.0)*fac(13.0))*Rm0n13(x, y, z); }

double RRm14n14(double x, double y, double z) { return sqrt(fac(0.0)*fac(28.0))*Rm14n14(x, y, z); }
double RRm_14n14(double x, double y, double z) { return sqrt(fac(28.0)*fac(0.0))*Rm_14n14(x, y, z); }
double RRm13n14(double x, double y, double z) { return sqrt(fac(1.0)*fac(27.0))*Rm13n14(x, y, z); }
double RRm_13n14(double x, double y, double z) { return sqrt(fac(27.0)*fac(1.0))*Rm_13n14(x, y, z); }
double RRm12n14(double x, double y, double z) { return sqrt(fac(2.0)*fac(26.0))*Rm12n14(x, y, z); }
double RRm_12n14(double x, double y, double z) { return sqrt(fac(26.0)*fac(2.0))*Rm_12n14(x, y, z); }
double RRm11n14(double x, double y, double z) { return sqrt(fac(3.0)*fac(25.0))*Rm11n14(x, y, z); }
double RRm_11n14(double x, double y, double z) { return sqrt(fac(25.0)*fac(3.0))*Rm_11n14(x, y, z); }
double RRm10n14(double x, double y, double z) { return sqrt(fac(4.0)*fac(24.0))*Rm10n14(x, y, z); }
double RRm_10n14(double x, double y, double z) { return sqrt(fac(24.0)*fac(4.0))*Rm_10n14(x, y, z); }
double RRm9n14(double x, double y, double z) { return sqrt(fac(5.0)*fac(23.0))*Rm9n14(x, y, z); }
double RRm_9n14(double x, double y, double z) { return sqrt(fac(23.0)*fac(5.0))*Rm_9n14(x, y, z); }
double RRm8n14(double x, double y, double z) { return sqrt(fac(6.0)*fac(22.0))*Rm8n14(x, y, z); }
double RRm_8n14(double x, double y, double z) { return sqrt(fac(22.0)*fac(6.0))*Rm_8n14(x, y, z); }
double RRm7n14(double x, double y, double z) { return sqrt(fac(7.0)*fac(21.0))*Rm7n14(x, y, z); }
double RRm_7n14(double x, double y, double z) { return sqrt(fac(21.0)*fac(7.0))*Rm_7n14(x, y, z); }
double RRm6n14(double x, double y, double z) { return sqrt(fac(8.0)*fac(20.0))*Rm6n14(x, y, z); }
double RRm_6n14(double x, double y, double z) { return sqrt(fac(20.0)*fac(8.0))*Rm_6n14(x, y, z); }
double RRm5n14(double x, double y, double z) { return sqrt(fac(9.0)*fac(19.0))*Rm5n14(x, y, z); }
double RRm_5n14(double x, double y, double z) { return sqrt(fac(19.0)*fac(9.0))*Rm_5n14(x, y, z); }
double RRm4n14(double x, double y, double z) { return sqrt(fac(10.0)*fac(18.0))*Rm4n14(x, y, z); }
double RRm_4n14(double x, double y, double z) { return sqrt(fac(18.0)*fac(10.0))*Rm_4n14(x, y, z); }
double RRm3n14(double x, double y, double z) { return sqrt(fac(11.0)*fac(17.0))*Rm3n14(x, y, z); }
double RRm_3n14(double x, double y, double z) { return sqrt(fac(17.0)*fac(11.0))*Rm_3n14(x, y, z); }
double RRm2n14(double x, double y, double z) { return sqrt(fac(12.0)*fac(16.0))*Rm2n14(x, y, z); }
double RRm_2n14(double x, double y, double z) { return sqrt(fac(16.0)*fac(12.0))*Rm_2n14(x, y, z); }
double RRm1n14(double x, double y, double z) { return sqrt(fac(13.0)*fac(15.0))*Rm1n14(x, y, z); }
double RRm_1n14(double x, double y, double z) { return sqrt(fac(15.0)*fac(13.0))*Rm_1n14(x, y, z); }
double RRm0n14(double x, double y, double z) { return sqrt(fac(14.0)*fac(14.0))*Rm0n14(x, y, z); }

double RRm15n15(double x, double y, double z) { return sqrt(fac(0.0)*fac(30.0))*Rm15n15(x, y, z); }
double RRm_15n15(double x, double y, double z) { return sqrt(fac(30.0)*fac(0.0))*Rm_15n15(x, y, z); }
double RRm14n15(double x, double y, double z) { return sqrt(fac(1.0)*fac(29.0))*Rm14n15(x, y, z); }
double RRm_14n15(double x, double y, double z) { return sqrt(fac(29.0)*fac(1.0))*Rm_14n15(x, y, z); }
double RRm13n15(double x, double y, double z) { return sqrt(fac(2.0)*fac(28.0))*Rm13n15(x, y, z); }
double RRm_13n15(double x, double y, double z) { return sqrt(fac(28.0)*fac(2.0))*Rm_13n15(x, y, z); }
double RRm12n15(double x, double y, double z) { return sqrt(fac(3.0)*fac(27.0))*Rm12n15(x, y, z); }
double RRm_12n15(double x, double y, double z) { return sqrt(fac(27.0)*fac(3.0))*Rm_12n15(x, y, z); }
double RRm11n15(double x, double y, double z) { return sqrt(fac(4.0)*fac(26.0))*Rm11n15(x, y, z); }
double RRm_11n15(double x, double y, double z) { return sqrt(fac(26.0)*fac(4.0))*Rm_11n15(x, y, z); }
double RRm10n15(double x, double y, double z) { return sqrt(fac(5.0)*fac(25.0))*Rm10n15(x, y, z); }
double RRm_10n15(double x, double y, double z) { return sqrt(fac(25.0)*fac(5.0))*Rm_10n15(x, y, z); }
double RRm9n15(double x, double y, double z) { return sqrt(fac(6.0)*fac(24.0))*Rm9n15(x, y, z); }
double RRm_9n15(double x, double y, double z) { return sqrt(fac(24.0)*fac(6.0))*Rm_9n15(x, y, z); }
double RRm8n15(double x, double y, double z) { return sqrt(fac(7.0)*fac(23.0))*Rm8n15(x, y, z); }
double RRm_8n15(double x, double y, double z) { return sqrt(fac(23.0)*fac(7.0))*Rm_8n15(x, y, z); }
double RRm7n15(double x, double y, double z) { return sqrt(fac(8.0)*fac(22.0))*Rm7n15(x, y, z); }
double RRm_7n15(double x, double y, double z) { return sqrt(fac(22.0)*fac(8.0))*Rm_7n15(x, y, z); }
double RRm6n15(double x, double y, double z) { return sqrt(fac(9.0)*fac(21.0))*Rm6n15(x, y, z); }
double RRm_6n15(double x, double y, double z) { return sqrt(fac(21.0)*fac(9.0))*Rm_6n15(x, y, z); }
double RRm5n15(double x, double y, double z) { return sqrt(fac(10.0)*fac(20.0))*Rm5n15(x, y, z); }
double RRm_5n15(double x, double y, double z) { return sqrt(fac(20.0)*fac(10.0))*Rm_5n15(x, y, z); }
double RRm4n15(double x, double y, double z) { return sqrt(fac(11.0)*fac(19.0))*Rm4n15(x, y, z); }
double RRm_4n15(double x, double y, double z) { return sqrt(fac(19.0)*fac(11.0))*Rm_4n15(x, y, z); }
double RRm3n15(double x, double y, double z) { return sqrt(fac(12.0)*fac(18.0))*Rm3n15(x, y, z); }
double RRm_3n15(double x, double y, double z) { return sqrt(fac(18.0)*fac(12.0))*Rm_3n15(x, y, z); }
double RRm2n15(double x, double y, double z) { return sqrt(fac(13.0)*fac(17.0))*Rm2n15(x, y, z); }
double RRm_2n15(double x, double y, double z) { return sqrt(fac(17.0)*fac(13.0))*Rm_2n15(x, y, z); }
double RRm1n15(double x, double y, double z) { return sqrt(fac(14.0)*fac(16.0))*Rm1n15(x, y, z); }
double RRm_1n15(double x, double y, double z) { return sqrt(fac(16.0)*fac(14.0))*Rm_1n15(x, y, z); }
double RRm0n15(double x, double y, double z) { return sqrt(fac(15.0)*fac(15.0))*Rm0n15(x, y, z); }

double RRm16n16(double x, double y, double z) { return sqrt(fac(0.0)*fac(32.0))*Rm16n16(x, y, z); }
double RRm_16n16(double x, double y, double z) { return sqrt(fac(32.0)*fac(0.0))*Rm_16n16(x, y, z); }
double RRm15n16(double x, double y, double z) { return sqrt(fac(1.0)*fac(31.0))*Rm15n16(x, y, z); }
double RRm_15n16(double x, double y, double z) { return sqrt(fac(31.0)*fac(1.0))*Rm_15n16(x, y, z); }
double RRm14n16(double x, double y, double z) { return sqrt(fac(2.0)*fac(30.0))*Rm14n16(x, y, z); }
double RRm_14n16(double x, double y, double z) { return sqrt(fac(30.0)*fac(2.0))*Rm_14n16(x, y, z); }
double RRm13n16(double x, double y, double z) { return sqrt(fac(3.0)*fac(29.0))*Rm13n16(x, y, z); }
double RRm_13n16(double x, double y, double z) { return sqrt(fac(29.0)*fac(3.0))*Rm_13n16(x, y, z); }
double RRm12n16(double x, double y, double z) { return sqrt(fac(4.0)*fac(28.0))*Rm12n16(x, y, z); }
double RRm_12n16(double x, double y, double z) { return sqrt(fac(28.0)*fac(4.0))*Rm_12n16(x, y, z); }
double RRm11n16(double x, double y, double z) { return sqrt(fac(5.0)*fac(27.0))*Rm11n16(x, y, z); }
double RRm_11n16(double x, double y, double z) { return sqrt(fac(27.0)*fac(5.0))*Rm_11n16(x, y, z); }
double RRm10n16(double x, double y, double z) { return sqrt(fac(6.0)*fac(26.0))*Rm10n16(x, y, z); }
double RRm_10n16(double x, double y, double z) { return sqrt(fac(26.0)*fac(6.0))*Rm_10n16(x, y, z); }
double RRm9n16(double x, double y, double z) { return sqrt(fac(7.0)*fac(25.0))*Rm9n16(x, y, z); }
double RRm_9n16(double x, double y, double z) { return sqrt(fac(25.0)*fac(7.0))*Rm_9n16(x, y, z); }
double RRm8n16(double x, double y, double z) { return sqrt(fac(8.0)*fac(24.0))*Rm8n16(x, y, z); }
double RRm_8n16(double x, double y, double z) { return sqrt(fac(24.0)*fac(8.0))*Rm_8n16(x, y, z); }
double RRm7n16(double x, double y, double z) { return sqrt(fac(9.0)*fac(23.0))*Rm7n16(x, y, z); }
double RRm_7n16(double x, double y, double z) { return sqrt(fac(23.0)*fac(9.0))*Rm_7n16(x, y, z); }
double RRm6n16(double x, double y, double z) { return sqrt(fac(10.0)*fac(22.0))*Rm6n16(x, y, z); }
double RRm_6n16(double x, double y, double z) { return sqrt(fac(22.0)*fac(10.0))*Rm_6n16(x, y, z); }
double RRm5n16(double x, double y, double z) { return sqrt(fac(11.0)*fac(21.0))*Rm5n16(x, y, z); }
double RRm_5n16(double x, double y, double z) { return sqrt(fac(21.0)*fac(11.0))*Rm_5n16(x, y, z); }
double RRm4n16(double x, double y, double z) { return sqrt(fac(12.0)*fac(20.0))*Rm4n16(x, y, z); }
double RRm_4n16(double x, double y, double z) { return sqrt(fac(20.0)*fac(12.0))*Rm_4n16(x, y, z); }
double RRm3n16(double x, double y, double z) { return sqrt(fac(13.0)*fac(19.0))*Rm3n16(x, y, z); }
double RRm_3n16(double x, double y, double z) { return sqrt(fac(19.0)*fac(13.0))*Rm_3n16(x, y, z); }
double RRm2n16(double x, double y, double z) { return sqrt(fac(14.0)*fac(18.0))*Rm2n16(x, y, z); }
double RRm_2n16(double x, double y, double z) { return sqrt(fac(18.0)*fac(14.0))*Rm_2n16(x, y, z); }
double RRm1n16(double x, double y, double z) { return sqrt(fac(15.0)*fac(17.0))*Rm1n16(x, y, z); }
double RRm_1n16(double x, double y, double z) { return sqrt(fac(17.0)*fac(15.0))*Rm_1n16(x, y, z); }
double RRm0n16(double x, double y, double z) { return sqrt(fac(16.0)*fac(16.0))*Rm0n16(x, y, z); }

double RRm17n17(double x, double y, double z) { return sqrt(fac(0.0)*fac(34.0))*Rm17n17(x, y, z); }
double RRm_17n17(double x, double y, double z) { return sqrt(fac(34.0)*fac(0.0))*Rm_17n17(x, y, z); }
double RRm16n17(double x, double y, double z) { return sqrt(fac(1.0)*fac(33.0))*Rm16n17(x, y, z); }
double RRm_16n17(double x, double y, double z) { return sqrt(fac(33.0)*fac(1.0))*Rm_16n17(x, y, z); }
double RRm15n17(double x, double y, double z) { return sqrt(fac(2.0)*fac(32.0))*Rm15n17(x, y, z); }
double RRm_15n17(double x, double y, double z) { return sqrt(fac(32.0)*fac(2.0))*Rm_15n17(x, y, z); }
double RRm14n17(double x, double y, double z) { return sqrt(fac(3.0)*fac(31.0))*Rm14n17(x, y, z); }
double RRm_14n17(double x, double y, double z) { return sqrt(fac(31.0)*fac(3.0))*Rm_14n17(x, y, z); }
double RRm13n17(double x, double y, double z) { return sqrt(fac(4.0)*fac(30.0))*Rm13n17(x, y, z); }
double RRm_13n17(double x, double y, double z) { return sqrt(fac(30.0)*fac(4.0))*Rm_13n17(x, y, z); }
double RRm12n17(double x, double y, double z) { return sqrt(fac(5.0)*fac(29.0))*Rm12n17(x, y, z); }
double RRm_12n17(double x, double y, double z) { return sqrt(fac(29.0)*fac(5.0))*Rm_12n17(x, y, z); }
double RRm11n17(double x, double y, double z) { return sqrt(fac(6.0)*fac(28.0))*Rm11n17(x, y, z); }
double RRm_11n17(double x, double y, double z) { return sqrt(fac(28.0)*fac(6.0))*Rm_11n17(x, y, z); }
double RRm10n17(double x, double y, double z) { return sqrt(fac(7.0)*fac(27.0))*Rm10n17(x, y, z); }
double RRm_10n17(double x, double y, double z) { return sqrt(fac(27.0)*fac(7.0))*Rm_10n17(x, y, z); }
double RRm9n17(double x, double y, double z) { return sqrt(fac(8.0)*fac(26.0))*Rm9n17(x, y, z); }
double RRm_9n17(double x, double y, double z) { return sqrt(fac(26.0)*fac(8.0))*Rm_9n17(x, y, z); }
double RRm8n17(double x, double y, double z) { return sqrt(fac(9.0)*fac(25.0))*Rm8n17(x, y, z); }
double RRm_8n17(double x, double y, double z) { return sqrt(fac(25.0)*fac(9.0))*Rm_8n17(x, y, z); }
double RRm7n17(double x, double y, double z) { return sqrt(fac(10.0)*fac(24.0))*Rm7n17(x, y, z); }
double RRm_7n17(double x, double y, double z) { return sqrt(fac(24.0)*fac(10.0))*Rm_7n17(x, y, z); }
double RRm6n17(double x, double y, double z) { return sqrt(fac(11.0)*fac(23.0))*Rm6n17(x, y, z); }
double RRm_6n17(double x, double y, double z) { return sqrt(fac(23.0)*fac(11.0))*Rm_6n17(x, y, z); }
double RRm5n17(double x, double y, double z) { return sqrt(fac(12.0)*fac(22.0))*Rm5n17(x, y, z); }
double RRm_5n17(double x, double y, double z) { return sqrt(fac(22.0)*fac(12.0))*Rm_5n17(x, y, z); }
double RRm4n17(double x, double y, double z) { return sqrt(fac(13.0)*fac(21.0))*Rm4n17(x, y, z); }
double RRm_4n17(double x, double y, double z) { return sqrt(fac(21.0)*fac(13.0))*Rm_4n17(x, y, z); }
double RRm3n17(double x, double y, double z) { return sqrt(fac(14.0)*fac(20.0))*Rm3n17(x, y, z); }
double RRm_3n17(double x, double y, double z) { return sqrt(fac(20.0)*fac(14.0))*Rm_3n17(x, y, z); }
double RRm2n17(double x, double y, double z) { return sqrt(fac(15.0)*fac(19.0))*Rm2n17(x, y, z); }
double RRm_2n17(double x, double y, double z) { return sqrt(fac(19.0)*fac(15.0))*Rm_2n17(x, y, z); }
double RRm1n17(double x, double y, double z) { return sqrt(fac(16.0)*fac(18.0))*Rm1n17(x, y, z); }
double RRm_1n17(double x, double y, double z) { return sqrt(fac(18.0)*fac(16.0))*Rm_1n17(x, y, z); }
double RRm0n17(double x, double y, double z) { return sqrt(fac(17.0)*fac(17.0))*Rm0n17(x, y, z); }

double RRm18n18(double x, double y, double z) { return sqrt(fac(0.0)*fac(36.0))*Rm18n18(x, y, z); }
double RRm_18n18(double x, double y, double z) { return sqrt(fac(36.0)*fac(0.0))*Rm_18n18(x, y, z); }
double RRm17n18(double x, double y, double z) { return sqrt(fac(1.0)*fac(35.0))*Rm17n18(x, y, z); }
double RRm_17n18(double x, double y, double z) { return sqrt(fac(35.0)*fac(1.0))*Rm_17n18(x, y, z); }
double RRm16n18(double x, double y, double z) { return sqrt(fac(2.0)*fac(34.0))*Rm16n18(x, y, z); }
double RRm_16n18(double x, double y, double z) { return sqrt(fac(34.0)*fac(2.0))*Rm_16n18(x, y, z); }
double RRm15n18(double x, double y, double z) { return sqrt(fac(3.0)*fac(33.0))*Rm15n18(x, y, z); }
double RRm_15n18(double x, double y, double z) { return sqrt(fac(33.0)*fac(3.0))*Rm_15n18(x, y, z); }
double RRm14n18(double x, double y, double z) { return sqrt(fac(4.0)*fac(32.0))*Rm14n18(x, y, z); }
double RRm_14n18(double x, double y, double z) { return sqrt(fac(32.0)*fac(4.0))*Rm_14n18(x, y, z); }
double RRm13n18(double x, double y, double z) { return sqrt(fac(5.0)*fac(31.0))*Rm13n18(x, y, z); }
double RRm_13n18(double x, double y, double z) { return sqrt(fac(31.0)*fac(5.0))*Rm_13n18(x, y, z); }
double RRm12n18(double x, double y, double z) { return sqrt(fac(6.0)*fac(30.0))*Rm12n18(x, y, z); }
double RRm_12n18(double x, double y, double z) { return sqrt(fac(30.0)*fac(6.0))*Rm_12n18(x, y, z); }
double RRm11n18(double x, double y, double z) { return sqrt(fac(7.0)*fac(29.0))*Rm11n18(x, y, z); }
double RRm_11n18(double x, double y, double z) { return sqrt(fac(29.0)*fac(7.0))*Rm_11n18(x, y, z); }
double RRm10n18(double x, double y, double z) { return sqrt(fac(8.0)*fac(28.0))*Rm10n18(x, y, z); }
double RRm_10n18(double x, double y, double z) { return sqrt(fac(28.0)*fac(8.0))*Rm_10n18(x, y, z); }
double RRm9n18(double x, double y, double z) { return sqrt(fac(9.0)*fac(27.0))*Rm9n18(x, y, z); }
double RRm_9n18(double x, double y, double z) { return sqrt(fac(27.0)*fac(9.0))*Rm_9n18(x, y, z); }
double RRm8n18(double x, double y, double z) { return sqrt(fac(10.0)*fac(26.0))*Rm8n18(x, y, z); }
double RRm_8n18(double x, double y, double z) { return sqrt(fac(26.0)*fac(10.0))*Rm_8n18(x, y, z); }
double RRm7n18(double x, double y, double z) { return sqrt(fac(11.0)*fac(25.0))*Rm7n18(x, y, z); }
double RRm_7n18(double x, double y, double z) { return sqrt(fac(25.0)*fac(11.0))*Rm_7n18(x, y, z); }
double RRm6n18(double x, double y, double z) { return sqrt(fac(12.0)*fac(24.0))*Rm6n18(x, y, z); }
double RRm_6n18(double x, double y, double z) { return sqrt(fac(24.0)*fac(12.0))*Rm_6n18(x, y, z); }
double RRm5n18(double x, double y, double z) { return sqrt(fac(13.0)*fac(23.0))*Rm5n18(x, y, z); }
double RRm_5n18(double x, double y, double z) { return sqrt(fac(23.0)*fac(13.0))*Rm_5n18(x, y, z); }
double RRm4n18(double x, double y, double z) { return sqrt(fac(14.0)*fac(22.0))*Rm4n18(x, y, z); }
double RRm_4n18(double x, double y, double z) { return sqrt(fac(22.0)*fac(14.0))*Rm_4n18(x, y, z); }
double RRm3n18(double x, double y, double z) { return sqrt(fac(15.0)*fac(21.0))*Rm3n18(x, y, z); }
double RRm_3n18(double x, double y, double z) { return sqrt(fac(21.0)*fac(15.0))*Rm_3n18(x, y, z); }
double RRm2n18(double x, double y, double z) { return sqrt(fac(16.0)*fac(20.0))*Rm2n18(x, y, z); }
double RRm_2n18(double x, double y, double z) { return sqrt(fac(20.0)*fac(16.0))*Rm_2n18(x, y, z); }
double RRm1n18(double x, double y, double z) { return sqrt(fac(17.0)*fac(19.0))*Rm1n18(x, y, z); }
double RRm_1n18(double x, double y, double z) { return sqrt(fac(19.0)*fac(17.0))*Rm_1n18(x, y, z); }
double RRm0n18(double x, double y, double z) { return sqrt(fac(18.0)*fac(18.0))*Rm0n18(x, y, z); }

double RRm19n19(double x, double y, double z) { return sqrt(fac(0.0)*fac(38.0))*Rm19n19(x, y, z); }
double RRm_19n19(double x, double y, double z) { return sqrt(fac(38.0)*fac(0.0))*Rm_19n19(x, y, z); }
double RRm18n19(double x, double y, double z) { return sqrt(fac(1.0)*fac(37.0))*Rm18n19(x, y, z); }
double RRm_18n19(double x, double y, double z) { return sqrt(fac(37.0)*fac(1.0))*Rm_18n19(x, y, z); }
double RRm17n19(double x, double y, double z) { return sqrt(fac(2.0)*fac(36.0))*Rm17n19(x, y, z); }
double RRm_17n19(double x, double y, double z) { return sqrt(fac(36.0)*fac(2.0))*Rm_17n19(x, y, z); }
double RRm16n19(double x, double y, double z) { return sqrt(fac(3.0)*fac(35.0))*Rm16n19(x, y, z); }
double RRm_16n19(double x, double y, double z) { return sqrt(fac(35.0)*fac(3.0))*Rm_16n19(x, y, z); }
double RRm15n19(double x, double y, double z) { return sqrt(fac(4.0)*fac(34.0))*Rm15n19(x, y, z); }
double RRm_15n19(double x, double y, double z) { return sqrt(fac(34.0)*fac(4.0))*Rm_15n19(x, y, z); }
double RRm14n19(double x, double y, double z) { return sqrt(fac(5.0)*fac(33.0))*Rm14n19(x, y, z); }
double RRm_14n19(double x, double y, double z) { return sqrt(fac(33.0)*fac(5.0))*Rm_14n19(x, y, z); }
double RRm13n19(double x, double y, double z) { return sqrt(fac(6.0)*fac(32.0))*Rm13n19(x, y, z); }
double RRm_13n19(double x, double y, double z) { return sqrt(fac(32.0)*fac(6.0))*Rm_13n19(x, y, z); }
double RRm12n19(double x, double y, double z) { return sqrt(fac(7.0)*fac(31.0))*Rm12n19(x, y, z); }
double RRm_12n19(double x, double y, double z) { return sqrt(fac(31.0)*fac(7.0))*Rm_12n19(x, y, z); }
double RRm11n19(double x, double y, double z) { return sqrt(fac(8.0)*fac(30.0))*Rm11n19(x, y, z); }
double RRm_11n19(double x, double y, double z) { return sqrt(fac(30.0)*fac(8.0))*Rm_11n19(x, y, z); }
double RRm10n19(double x, double y, double z) { return sqrt(fac(9.0)*fac(29.0))*Rm10n19(x, y, z); }
double RRm_10n19(double x, double y, double z) { return sqrt(fac(29.0)*fac(9.0))*Rm_10n19(x, y, z); }
double RRm9n19(double x, double y, double z) { return sqrt(fac(10.0)*fac(28.0))*Rm9n19(x, y, z); }
double RRm_9n19(double x, double y, double z) { return sqrt(fac(28.0)*fac(10.0))*Rm_9n19(x, y, z); }
double RRm8n19(double x, double y, double z) { return sqrt(fac(11.0)*fac(27.0))*Rm8n19(x, y, z); }
double RRm_8n19(double x, double y, double z) { return sqrt(fac(27.0)*fac(11.0))*Rm_8n19(x, y, z); }
double RRm7n19(double x, double y, double z) { return sqrt(fac(12.0)*fac(26.0))*Rm7n19(x, y, z); }
double RRm_7n19(double x, double y, double z) { return sqrt(fac(26.0)*fac(12.0))*Rm_7n19(x, y, z); }
double RRm6n19(double x, double y, double z) { return sqrt(fac(13.0)*fac(25.0))*Rm6n19(x, y, z); }
double RRm_6n19(double x, double y, double z) { return sqrt(fac(25.0)*fac(13.0))*Rm_6n19(x, y, z); }
double RRm5n19(double x, double y, double z) { return sqrt(fac(14.0)*fac(24.0))*Rm5n19(x, y, z); }
double RRm_5n19(double x, double y, double z) { return sqrt(fac(24.0)*fac(14.0))*Rm_5n19(x, y, z); }
double RRm4n19(double x, double y, double z) { return sqrt(fac(15.0)*fac(23.0))*Rm4n19(x, y, z); }
double RRm_4n19(double x, double y, double z) { return sqrt(fac(23.0)*fac(15.0))*Rm_4n19(x, y, z); }
double RRm3n19(double x, double y, double z) { return sqrt(fac(16.0)*fac(22.0))*Rm3n19(x, y, z); }
double RRm_3n19(double x, double y, double z) { return sqrt(fac(22.0)*fac(16.0))*Rm_3n19(x, y, z); }
double RRm2n19(double x, double y, double z) { return sqrt(fac(17.0)*fac(21.0))*Rm2n19(x, y, z); }
double RRm_2n19(double x, double y, double z) { return sqrt(fac(21.0)*fac(17.0))*Rm_2n19(x, y, z); }
double RRm1n19(double x, double y, double z) { return sqrt(fac(18.0)*fac(20.0))*Rm1n19(x, y, z); }
double RRm_1n19(double x, double y, double z) { return sqrt(fac(20.0)*fac(18.0))*Rm_1n19(x, y, z); }
double RRm0n19(double x, double y, double z) { return sqrt(fac(19.0)*fac(19.0))*Rm0n19(x, y, z); }

//double (*P[])(double,double,double) = {Rm0n0,Rm_1n1,Rm0n1,Rm1n1,Rm_2n2,Rm_1n2,Rm0n2,Rm1n2,Rm2n2,Rm_3n3,Rm_2n3,Rm_1n3,Rm0n3,Rm1n3,Rm2n3,Rm3n3 };//Bases without normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=3
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=4
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3 };// Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=5
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=6
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=7
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=8
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=9
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=10
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2
//double(*P[])(double, double, double) = { Rm0n0,Rm_1n1,Rm0n1,Rm1n1,Rm_2n2,Rm_1n2,Rm0n2,Rm1n2,Rm2n2,Rm_3n3,Rm_2n3,Rm_1n3,Rm0n3,Rm1n3,Rm2n3,Rm3n3,Rm_4n4,Rm_3n4,Rm_2n4,Rm_1n4,Rm0n4,Rm1n4,Rm2n4,Rm3n4,Rm4n4,Rm_5n5,Rm_4n5,Rm_3n5,Rm_2n5,Rm_1n5,Rm0n5,Rm1n5,Rm2n5,Rm3n5,Rm4n5,Rm5n5 ,Rm_6n6,Rm_5n6,Rm_4n6,Rm_3n6,Rm_2n6,Rm_1n6,Rm0n6,Rm1n6,Rm2n6,Rm3n6,Rm4n6,Rm5n6,Rm6n6,Rm_7n7,Rm_6n7,Rm_5n7,Rm_4n7,Rm_3n7,Rm_2n7,Rm_1n7,Rm0n7,Rm1n7,Rm2n7,Rm3n7,Rm4n7,Rm5n7,Rm6n7,Rm7n7,Rm_8n8,Rm_7n8,Rm_6n8,Rm_5n8,Rm_4n8,Rm_3n8,Rm_2n8,Rm_1n8,Rm0n8,Rm1n8,Rm2n8,Rm3n8,Rm4n8,Rm5n8,Rm6n8,Rm7n8,Rm8n8,Rm_9n9,Rm_8n9,Rm_7n9,Rm_6n9,Rm_5n9,Rm_4n9,Rm_3n9,Rm_2n9,Rm_1n9,Rm0n9,Rm1n9,Rm2n9,Rm3n9,Rm4n9,Rm5n9,Rm6n9,Rm7n9,Rm8n9,Rm9n9 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=11
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10};//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=12
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=13
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=14
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=15
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=16
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14,RRm_15n15,RRm_14n15,RRm_13n15, RRm_12n15, RRm_11n15, RRm_10n15, RRm_9n15, RRm_8n15, RRm_7n15, RRm_6n15, RRm_5n15, RRm_4n15, RRm_3n15, RRm_2n15,RRm_1n15,RRm0n15,RRm1n15,RRm2n15,RRm3n15,RRm4n15,RRm5n15,RRm6n15,RRm7n15,RRm8n15,RRm9n15,RRm10n15,RRm11n15,RRm12n15,RRm13n15,RRm14n15,RRm15n15 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=17
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14,RRm_15n15,RRm_14n15,RRm_13n15, RRm_12n15, RRm_11n15, RRm_10n15, RRm_9n15, RRm_8n15, RRm_7n15, RRm_6n15, RRm_5n15, RRm_4n15, RRm_3n15, RRm_2n15,RRm_1n15,RRm0n15,RRm1n15,RRm2n15,RRm3n15,RRm4n15,RRm5n15,RRm6n15,RRm7n15,RRm8n15,RRm9n15,RRm10n15,RRm11n15,RRm12n15,RRm13n15,RRm14n15,RRm15n15,RRm_16n16,RRm_15n16,RRm_14n16,RRm_13n16, RRm_12n16, RRm_11n16, RRm_10n16, RRm_9n16, RRm_8n16, RRm_7n16, RRm_6n16, RRm_5n16, RRm_4n16, RRm_3n16, RRm_2n16,RRm_1n16,RRm0n16,RRm1n16,RRm2n16,RRm3n16,RRm4n16,RRm5n16,RRm6n16,RRm7n16,RRm8n16,RRm9n16,RRm10n16,RRm11n16,RRm12n16,RRm13n16,RRm14n16,RRm15n16,RRm16n16 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=18
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14,RRm_15n15,RRm_14n15,RRm_13n15, RRm_12n15, RRm_11n15, RRm_10n15, RRm_9n15, RRm_8n15, RRm_7n15, RRm_6n15, RRm_5n15, RRm_4n15, RRm_3n15, RRm_2n15,RRm_1n15,RRm0n15,RRm1n15,RRm2n15,RRm3n15,RRm4n15,RRm5n15,RRm6n15,RRm7n15,RRm8n15,RRm9n15,RRm10n15,RRm11n15,RRm12n15,RRm13n15,RRm14n15,RRm15n15,RRm_16n16,RRm_15n16,RRm_14n16,RRm_13n16, RRm_12n16, RRm_11n16, RRm_10n16, RRm_9n16, RRm_8n16, RRm_7n16, RRm_6n16, RRm_5n16, RRm_4n16, RRm_3n16, RRm_2n16,RRm_1n16,RRm0n16,RRm1n16,RRm2n16,RRm3n16,RRm4n16,RRm5n16,RRm6n16,RRm7n16,RRm8n16,RRm9n16,RRm10n16,RRm11n16,RRm12n16,RRm13n16,RRm14n16,RRm15n16,RRm16n16,RRm_17n17,RRm_16n17,RRm_15n17,RRm_14n17,RRm_13n17, RRm_12n17, RRm_11n17, RRm_10n17, RRm_9n17, RRm_8n17, RRm_7n17, RRm_6n17, RRm_5n17, RRm_4n17, RRm_3n17, RRm_2n17,RRm_1n17,RRm0n17,RRm1n17,RRm2n17,RRm3n17,RRm4n17,RRm5n17,RRm6n17,RRm7n17,RRm8n17,RRm9n17,RRm10n17,RRm11n17,RRm12n17,RRm13n17,RRm14n17,RRm15n17,RRm16n17,RRm17n17 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=19
//double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14,RRm_15n15,RRm_14n15,RRm_13n15, RRm_12n15, RRm_11n15, RRm_10n15, RRm_9n15, RRm_8n15, RRm_7n15, RRm_6n15, RRm_5n15, RRm_4n15, RRm_3n15, RRm_2n15,RRm_1n15,RRm0n15,RRm1n15,RRm2n15,RRm3n15,RRm4n15,RRm5n15,RRm6n15,RRm7n15,RRm8n15,RRm9n15,RRm10n15,RRm11n15,RRm12n15,RRm13n15,RRm14n15,RRm15n15,RRm_16n16,RRm_15n16,RRm_14n16,RRm_13n16, RRm_12n16, RRm_11n16, RRm_10n16, RRm_9n16, RRm_8n16, RRm_7n16, RRm_6n16, RRm_5n16, RRm_4n16, RRm_3n16, RRm_2n16,RRm_1n16,RRm0n16,RRm1n16,RRm2n16,RRm3n16,RRm4n16,RRm5n16,RRm6n16,RRm7n16,RRm8n16,RRm9n16,RRm10n16,RRm11n16,RRm12n16,RRm13n16,RRm14n16,RRm15n16,RRm16n16,RRm_17n17,RRm_16n17,RRm_15n17,RRm_14n17,RRm_13n17, RRm_12n17, RRm_11n17, RRm_10n17, RRm_9n17, RRm_8n17, RRm_7n17, RRm_6n17, RRm_5n17, RRm_4n17, RRm_3n17, RRm_2n17,RRm_1n17,RRm0n17,RRm1n17,RRm2n17,RRm3n17,RRm4n17,RRm5n17,RRm6n17,RRm7n17,RRm8n17,RRm9n17,RRm10n17,RRm11n17,RRm12n17,RRm13n17,RRm14n17,RRm15n17,RRm16n17,RRm17n17,RRm_18n18,RRm_17n18,RRm_16n18,RRm_15n18,RRm_14n18,RRm_13n18, RRm_12n18, RRm_11n18, RRm_10n18, RRm_9n18, RRm_8n18, RRm_7n18, RRm_6n18, RRm_5n18, RRm_4n18, RRm_3n18, RRm_2n18,RRm_1n18,RRm0n18,RRm1n18,RRm2n18,RRm3n18,RRm4n18,RRm5n18,RRm6n18,RRm7n18,RRm8n18,RRm9n18,RRm10n18,RRm11n18,RRm12n18,RRm13n18,RRm14n18,RRm15n18,RRm16n18,RRm17n18,RRm18n18 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2

//p=20
double(*P[])(double, double, double) = { RRm0n0,RRm_1n1,RRm0n1,RRm1n1,RRm_2n2,RRm_1n2,RRm0n2,RRm1n2,RRm2n2,RRm_3n3,RRm_2n3,RRm_1n3,RRm0n3,RRm1n3,RRm2n3,RRm3n3,RRm_4n4,RRm_3n4,RRm_2n4,RRm_1n4,RRm0n4,RRm1n4,RRm2n4,RRm3n4,RRm4n4,RRm_5n5,RRm_4n5,RRm_3n5,RRm_2n5,RRm_1n5,RRm0n5,RRm1n5,RRm2n5,RRm3n5,RRm4n5,RRm5n5 ,RRm_6n6,RRm_5n6,RRm_4n6,RRm_3n6,RRm_2n6,RRm_1n6,RRm0n6,RRm1n6,RRm2n6,RRm3n6,RRm4n6,RRm5n6,RRm6n6,RRm_7n7,RRm_6n7,RRm_5n7,RRm_4n7,RRm_3n7,RRm_2n7,RRm_1n7,RRm0n7,RRm1n7,RRm2n7,RRm3n7,RRm4n7,RRm5n7,RRm6n7,RRm7n7,RRm_8n8,RRm_7n8,RRm_6n8,RRm_5n8,RRm_4n8,RRm_3n8,RRm_2n8,RRm_1n8,RRm0n8,RRm1n8,RRm2n8,RRm3n8,RRm4n8,RRm5n8,RRm6n8,RRm7n8,RRm8n8,RRm_9n9,RRm_8n9,RRm_7n9,RRm_6n9,RRm_5n9,RRm_4n9,RRm_3n9,RRm_2n9,RRm_1n9,RRm0n9,RRm1n9,RRm2n9,RRm3n9,RRm4n9,RRm5n9,RRm6n9,RRm7n9,RRm8n9,RRm9n9,RRm_10n10,RRm_9n10,RRm_8n10,RRm_7n10,RRm_6n10,RRm_5n10,RRm_4n10,RRm_3n10,RRm_2n10,RRm_1n10,RRm0n10,RRm1n10,RRm2n10,RRm3n10,RRm4n10,RRm5n10,RRm6n10,RRm7n10,RRm8n10,RRm9n10,RRm10n10,RRm_11n11,RRm_10n11,RRm_9n11,RRm_8n11,RRm_7n11,RRm_6n11,RRm_5n11,RRm_4n11,RRm_3n11,RRm_2n11,RRm_1n11,RRm0n11,RRm1n11,RRm2n11,RRm3n11,RRm4n11,RRm5n11,RRm6n11,RRm7n11,RRm8n11,RRm9n11,RRm10n11,RRm11n11,RRm_12n12, RRm_11n12, RRm_10n12, RRm_9n12, RRm_8n12, RRm_7n12, RRm_6n12, RRm_5n12, RRm_4n12, RRm_3n12, RRm_2n12, RRm_1n12,RRm0n12,RRm1n12,RRm2n12,RRm3n12,RRm4n12,RRm5n12,RRm6n12,RRm7n12,RRm8n12,RRm9n12,RRm10n12,RRm11n12,RRm12n12,RRm_13n13, RRm_12n13, RRm_11n13, RRm_10n13, RRm_9n13, RRm_8n13, RRm_7n13, RRm_6n13, RRm_5n13, RRm_4n13, RRm_3n13, RRm_2n13,RRm_1n13,RRm0n13,RRm1n13,RRm2n13,RRm3n13,RRm4n13,RRm5n13,RRm6n13,RRm7n13,RRm8n13,RRm9n13,RRm10n13,RRm11n13,RRm12n13,RRm13n13,RRm_14n14,RRm_13n14, RRm_12n14, RRm_11n14, RRm_10n14, RRm_9n14, RRm_8n14, RRm_7n14, RRm_6n14, RRm_5n14, RRm_4n14, RRm_3n14, RRm_2n14,RRm_1n14,RRm0n14,RRm1n14,RRm2n14,RRm3n14,RRm4n14,RRm5n14,RRm6n14,RRm7n14,RRm8n14,RRm9n14,RRm10n14,RRm11n14,RRm12n14,RRm13n14,RRm14n14,RRm_15n15,RRm_14n15,RRm_13n15, RRm_12n15, RRm_11n15, RRm_10n15, RRm_9n15, RRm_8n15, RRm_7n15, RRm_6n15, RRm_5n15, RRm_4n15, RRm_3n15, RRm_2n15,RRm_1n15,RRm0n15,RRm1n15,RRm2n15,RRm3n15,RRm4n15,RRm5n15,RRm6n15,RRm7n15,RRm8n15,RRm9n15,RRm10n15,RRm11n15,RRm12n15,RRm13n15,RRm14n15,RRm15n15,RRm_16n16,RRm_15n16,RRm_14n16,RRm_13n16, RRm_12n16, RRm_11n16, RRm_10n16, RRm_9n16, RRm_8n16, RRm_7n16, RRm_6n16, RRm_5n16, RRm_4n16, RRm_3n16, RRm_2n16,RRm_1n16,RRm0n16,RRm1n16,RRm2n16,RRm3n16,RRm4n16,RRm5n16,RRm6n16,RRm7n16,RRm8n16,RRm9n16,RRm10n16,RRm11n16,RRm12n16,RRm13n16,RRm14n16,RRm15n16,RRm16n16,RRm_17n17,RRm_16n17,RRm_15n17,RRm_14n17,RRm_13n17, RRm_12n17, RRm_11n17, RRm_10n17, RRm_9n17, RRm_8n17, RRm_7n17, RRm_6n17, RRm_5n17, RRm_4n17, RRm_3n17, RRm_2n17,RRm_1n17,RRm0n17,RRm1n17,RRm2n17,RRm3n17,RRm4n17,RRm5n17,RRm6n17,RRm7n17,RRm8n17,RRm9n17,RRm10n17,RRm11n17,RRm12n17,RRm13n17,RRm14n17,RRm15n17,RRm16n17,RRm17n17,RRm_18n18,RRm_17n18,RRm_16n18,RRm_15n18,RRm_14n18,RRm_13n18, RRm_12n18, RRm_11n18, RRm_10n18, RRm_9n18, RRm_8n18, RRm_7n18, RRm_6n18, RRm_5n18, RRm_4n18, RRm_3n18, RRm_2n18,RRm_1n18,RRm0n18,RRm1n18,RRm2n18,RRm3n18,RRm4n18,RRm5n18,RRm6n18,RRm7n18,RRm8n18,RRm9n18,RRm10n18,RRm11n18,RRm12n18,RRm13n18,RRm14n18,RRm15n18,RRm16n18,RRm17n18,RRm18n18,RRm_19n19,RRm_18n19,RRm_17n19,RRm_16n19,RRm_15n19,RRm_14n19,RRm_13n19, RRm_12n19, RRm_11n19, RRm_10n19, RRm_9n19, RRm_8n19, RRm_7n19, RRm_6n19, RRm_5n19, RRm_4n19, RRm_3n19, RRm_2n19,RRm_1n19,RRm0n19,RRm1n19,RRm2n19,RRm3n19,RRm4n19,RRm5n19,RRm6n19,RRm7n19,RRm8n19,RRm9n19,RRm10n19,RRm11n19,RRm12n19,RRm13n19,RRm14n19,RRm15n19,RRm16n19,RRm17n19,RRm18n19,RRm19n19 };//Bases with normalization. The number of bases are depending on small p. Big P =(small p)^2
*/