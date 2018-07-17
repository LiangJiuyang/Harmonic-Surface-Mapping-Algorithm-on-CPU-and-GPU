#include<iostream>
#include<fstream>
#include<math.h>
#include <cstdlib>
#include <cmath>
//#include<conio.h>
#include <limits>

#define pi 3.1415926535898
 
//Set size of center box
#define L1 2.0//length
#define L2 2.0//width
#define L3 2.0//height

//Set radius of sphere
#define R  (sqrt(L1*L1 + L2*L2 + L3*L3)) / 2.0

//Set quiz point
#define QuizNumber 400//number of quiz points
class QuizPointSum
{
public:
	double *x=new double[QuizNumber]; 
	double *y=new double[QuizNumber];
	double *z=new double[QuizNumber];
         //double x[QuizNumber];
                //double y[QuizNumber];
         //double z[QuizNumber];
	QuizPointSum()
	{
		for (int i = 0; i < QuizNumber; i++)
		{
			z[i] = (2 * (i + 1) - 1) / (QuizNumber + 0.00) - 1;
			x[i] = (sqrt(1 - z[i] * z[i])*cos(2 * pi*(i + 1)*0.618))*R;
			y[i] = (sqrt(1 - z[i] * z[i])*sin(2 * pi*(i + 1)*0.618))*R;
			z[i] = z[i] * R;
		}
	}
	QuizPointSum &operator=(const QuizPointSum &P);//operator overloading(exchange position only)
}Quiz;

QuizPointSum & QuizPointSum::operator=(const QuizPointSum & P)
{
	for (int i = 0; i < QuizNumber; i++)
	{
		this->x[i] = P.x[i];
		this->y[i] = P.y[i];
		this->z[i] = P.z[i];
	}
	return *this;
}

//Find corresponding points set of quiz points set in center box

QuizPointSum FindCorrespondingPoints(QuizPointSum Q)
{
	QuizPointSum T;//T and Q has the same constructor. So that if particle in Q is in center box, the same as T that we do not need to move it.  
	for (int i = 0; i < QuizNumber; i++)
	{
		if (abs(Q.x[i]) >= (L1 / 2))
		{
			T.x[i] = (abs(Q.x[i]) - L1*int((abs(Q.x[i]) + L1 / 2) / L1))*Q.x[i] / abs(Q.x[i]);
		}
		if (abs(Q.y[i]) >= (L2 / 2))
		{
			T.y[i] = (abs(Q.y[i]) - L2*int((abs(Q.y[i]) + L2 / 2) / L2))*Q.y[i] / abs(Q.y[i]);
		}
		if (abs(Q.z[i]) >= (L3 / 2))
		{
			T.z[i] = (abs(Q.z[i]) - L3*int((abs(Q.z[i]) + L3 / 2) / L3))*Q.z[i] / abs(Q.z[i]);
		}
	}
	return T;
}


