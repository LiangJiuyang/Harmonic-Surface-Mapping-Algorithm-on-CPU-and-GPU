#include"BuildChargeClass.h"
#include"BuildTestPoint.h"
#include<fstream>
#include<string>
using namespace std;
//#define N 64//Number of field source
//#define N 5000//Number of field source
#define CenterBoxSourceNumber 80*80*80 
   
/*
#define NSmallestBox 60
#define N NSmallestBox*NSmallestBox*NSmallestBox//Number of field source
#define times ((NSmallestBox+0.0))*((NSmallestBox+0.0))*((NSmallestBox+0.0))*((NSmallestBox+0.0))/4.0//multipole which the theoretical madelung constant should time
*/ 
    
//Set charge of field source  
//Particle *FieldCharge = new Particle[N];

 
void SetFieldCharge(Particle *FieldCharge)
{
         const int N=CenterBoxSourceNumber; 
	ifstream infile; 
	infile.open("q0_80_1.txt");      
	if (!infile) cout << "error" << endl;

	string str;  
	double t1; 

	//存入数组
	//cout << "存入数组" << endl;

	double(*a)[4] = (double(*)[4])malloc(N * sizeof(double[4]));
	double *p = &a[0][0];
	while (infile >> t1)             //遇到空白符结束
	{
		*p = t1;
		p++;
	}
	infile.close();

	for (int i = 0; i < N; i++)
	{
		FieldCharge[i].x = a[i][0];
		FieldCharge[i].y = a[i][1];
		FieldCharge[i].z = a[i][2];
		FieldCharge[i].q = a[i][3];
	}
	//for (int i = 0; i < 50; i++)
	//{
		//printf("%lf     %lf     %lf      %lf\n",FieldCharge[i].x,FieldCharge[i].y,FieldCharge[i].z,FieldCharge[i].q);
	//}
}

/*    Set for Madelung constant     */
/*
void SetFieldCharge(Particle *FieldCharge)
{
FieldCharge[0].x = 0.5;
FieldCharge[0].y = 0.5;
FieldCharge[0].z = 0.5;
FieldCharge[0].q = 1.0;

FieldCharge[1].x = 0.5;
FieldCharge[1].y = 0.5;
FieldCharge[1].z = -0.5;
FieldCharge[1].q = -1.0;

FieldCharge[2].x = 0.5;
FieldCharge[2].y = -0.5;
FieldCharge[2].z = 0.5;
FieldCharge[2].q = -1.0;

FieldCharge[3].x = -0.5;
FieldCharge[3].y = 0.5;
FieldCharge[3].z = 0.5;
FieldCharge[3].q = -1.0;

FieldCharge[4].x = 0.5;
FieldCharge[4].y = -0.5;
FieldCharge[4].z = -0.5;
FieldCharge[4].q = 1.0;

FieldCharge[5].x = -0.5;
FieldCharge[5].y = 0.5;
FieldCharge[5].z = -0.5;
FieldCharge[5].q = 1.0;

FieldCharge[6].x = -0.5;
FieldCharge[6].y = -0.5;
FieldCharge[6].z = 0.5;
FieldCharge[6].q = 1.0;

FieldCharge[7].x = -0.5;
FieldCharge[7].y = -0.5;
FieldCharge[7].z = -0.5;
FieldCharge[7].q = -1.0;
}
*/
/*
void SetFieldCharge(Particle *FieldCharge)
{
FieldCharge[0].x = 0.25;
FieldCharge[0].y = 0.25;
FieldCharge[0].z = 0.25;
FieldCharge[0].q = 1.0;

FieldCharge[1].x = 0.25;
FieldCharge[1].y = 0.25;
FieldCharge[1].z = 0.75;
FieldCharge[1].q = -1.0;

FieldCharge[2].x = 0.75;
FieldCharge[2].y = 0.25;
FieldCharge[2].z = 0.25;
FieldCharge[2].q = -1.0;

FieldCharge[3].x = 0.25;
FieldCharge[3].y = 0.75;
FieldCharge[3].z = 0.25;
FieldCharge[3].q = -1.0;

FieldCharge[4].x = 0.25;
FieldCharge[4].y = 0.75;
FieldCharge[4].z = 0.75;
FieldCharge[4].q = 1.0;

FieldCharge[5].x = 0.75;
FieldCharge[5].y = 0.25;
FieldCharge[5].z = 0.75;
FieldCharge[5].q = 1.0;

FieldCharge[6].x = 0.75;
FieldCharge[6].y = 0.75;
FieldCharge[6].z = 0.25;
FieldCharge[6].q = 1.0;

FieldCharge[7].x = 0.75;
FieldCharge[7].y = 0.75;
FieldCharge[7].z = 0.75;
FieldCharge[7].q = -1.0;

for (int j = 1; j < 8; j++)
{
for (int i = 0; i < 8; i++)
{
if (j == 1)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x - 1;
FieldCharge[j * 8 + i].y = FieldCharge[i].y - 1;
FieldCharge[j * 8 + i].z = FieldCharge[i].z - 1;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 2)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x - 1;
FieldCharge[j * 8 + i].y = FieldCharge[i].y ;
FieldCharge[j * 8 + i].z = FieldCharge[i].z - 1;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 3)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x ;
FieldCharge[j * 8 + i].y = FieldCharge[i].y - 1;
FieldCharge[j * 8 + i].z = FieldCharge[i].z - 1;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 4)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x - 1;
FieldCharge[j * 8 + i].y = FieldCharge[i].y - 1;
FieldCharge[j * 8 + i].z = FieldCharge[i].z ;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 5)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x ;
FieldCharge[j * 8 + i].y = FieldCharge[i].y ;
FieldCharge[j * 8 + i].z = FieldCharge[i].z - 1;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 6)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x - 1;
FieldCharge[j * 8 + i].y = FieldCharge[i].y ;
FieldCharge[j * 8 + i].z = FieldCharge[i].z ;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
if (j == 7)
{
FieldCharge[j * 8 + i].x = FieldCharge[i].x ;
FieldCharge[j * 8 + i].y = FieldCharge[i].y - 1 ;
FieldCharge[j * 8 + i].z = FieldCharge[i].z ;
FieldCharge[j * 8 + i].q = FieldCharge[i].q;
}
}
}
}
*/
/*
void SetFieldCharge(Particle *FieldCharge)
{
for (int i = 0; i < N; i++)
{
if (i < 2500)
{
FieldCharge[i].x = -0.9996 + i*0.0004;
FieldCharge[i].y = -0.9996 + i*0.0004;
FieldCharge[i].z = -0.9996 + i*0.0004;
}
if (i >= 2500)
{
FieldCharge[i].x = 0.0004 + (i - 2500)*0.0004;
FieldCharge[i].y = 0.0004 + (i - 2500)*0.0004;
FieldCharge[i].z = 0.0004 + (i - 2500)*0.0004;
}
if (i % 2 == 0)
{
FieldCharge[i].q = 1.0;
}
if(i % 2 == 1)
{
FieldCharge[i].q = -1.0;
}
}
}
*/


double fac(double t)//calculate factorial
{
	double s;
	if (abs(t - 1)<0.001 || abs(t)<0.001)
		s = 1.0;
	else
	{
		s = t*fac(t - 1) + 0.00;
	}
	return s;
}