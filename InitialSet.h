#include"BuildChargeClass.h"
#include"BuildTestPoint.h"
#include<fstream>
#include<string>
using namespace std;

//Number of field source
//#define N 64 //For Madelung constant, using this code
#define CenterBoxSourceNumber 80*80*80 //For large scale test, using this code
       
//Set charge of field source (For large scale example) 
void SetFieldCharge(Particle *FieldCharge)
{
         const int N=CenterBoxSourceNumber; 
	ifstream infile; 
	infile.open("q0_80_1.txt");      
	if (!infile) cout << "error" << endl;

	string str;  
	double t1; 

	double(*a)[4] = (double(*)[4])malloc(N * sizeof(double[4]));
	double *p = &a[0][0];
	while (infile >> t1)             //After each white space in the file, the ">>" operator will stop reading in the contents, until another >> operator is encountered.
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

//Calculate factorial
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
