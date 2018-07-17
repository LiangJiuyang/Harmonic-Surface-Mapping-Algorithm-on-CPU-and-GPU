#include<iostream>
#include<fstream>
#include<math.h>
#include <cstdlib>
#include <cmath>
//#include<conio.h>
#include <limits>
#include<iomanip>
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
using namespace std;
class Particle
{
public:
	double x,y,z;//position of particle
	double q;//electric charge of particle
	Particle &operator=(const Particle &P);//operator overloading(exchange position only)
}Par;

Particle &Particle::operator=(const Particle &P)
{
	this->x = P.x;
	this->y = P.y;
	this->z = P.z;
	this->q = P.q;
	return *this;
}