#include"CalculateNearFieldPotential.h"
double F = 144.0;//Number of projection points on Omega_b(with Radius Rb set in "设定外层球并给出所有场源点的坐标.h")
double Fp = 89.0;
int Np = F * 2.0 + 2.0;
double zz,zt,zs,phia,phib;
double deltaz = 2.0 / F;
double DeltaR = (R+0.0)/pow(10.0,6.0);//d6elta R==(abs(R'-R''))

class ProjectPointSum
{
public:
	double *x = new double[Np];
	double *y = new double[Np];
	double *z = new double[Np];
	double *q = new double[Np];
	ProjectPointSum()
	{
		for (int i = 0; i <= F; i++)
		{
		    zz = -1.0 + i * deltaz;
		    zt = zz + sin(pi * zz) / pi;
            zs = sqrt(1 - zt * zt);
            phia = pi * i * Fp / F;
            phib = pi + phia;
            //[,sin(phia) *zs,zt];
			x[2 * i] = cos(phia) * zs * (Rb+0.00);
			x[2 * i + 1] = cos(phib) * zs * (Rb+0.00);
			y[2 * i] = sin(phia) * zs * (Rb+0.00);
			y[2 * i + 1] = sin(phib)* zs * (Rb+0.00);
			z[2 * i] = zt * (Rb+0.00);
			z[2 * i + 1] = zt * (Rb+0.00);
			q[2 * i] = pi * deltaz * (1.0 + cos(pi * zz));
			q[2 * i + 1] = pi * deltaz * (1.0 + cos(pi * zz));
		}
	}
	ProjectPointSum &operator=(const ProjectPointSum &P);//operator overloading(exchange position only)
}Project;

//Thompson's distribution on the surface of sphere
/*
class ProjectPointSum
{
public:
	double *x=new double[Np];
	double *y=new double[Np];
	double *z=new double[Np];
	//double *q = new double[Np];
	ProjectPointSum()
	{
		for (int i = 0; i < Np; i++)
		{
			z[i] = (2 * (i + 1) - 1) / (Np + 0.00) - 1;
			x[i] = (sqrt(1 - z[i] * z[i])*cos(2 * pi*(i + 1)*0.618))*(Rb+0.0);
			y[i] = (sqrt(1 - z[i] * z[i])*sin(2 * pi*(i + 1)*0.618))*(Rb+0.0);
			z[i] = z[i] * (Rb+0.0);
		}
	}
	ProjectPointSum &operator=(const ProjectPointSum &P);//operator overloading(exchange position only)
}Project;
*/

ProjectPointSum & ProjectPointSum::operator=(const ProjectPointSum & P)
{
	for (int i = 0; i < Np; i++)
	{
		this->x[i] = P.x[i];
		this->y[i] = P.y[i];
		this->z[i] = P.z[i];
	}
	return *this;
}

double Distance(double x, double y, double z, double xx, double yy, double zz)//calculate distance
{
	double distance;
	distance = sqrt((x - xx)*(x - xx) + (y - yy)*(y - yy) + (z - zz)*(z - zz));
	return distance;
}

//This code is using multipole expansion and project
double CalculateFarFieldPotential(double *C, double x, double y, double z,ProjectPointSum D,ProjectPointSum DD)
{
	double Potential = 0.00;
	double Sigma_r;//Sigma_r is for sigma(r')
	double Sigma_rr;//Sigma_rr is for sigma(r'')
	double kkk;

	for (int i = 0; i < Np; i++)
	{
		double *Q=new double[p*p];
		CalculateMultipleExpansion(Q, D.x[i], D.y[i], D.z[i]);
		Sigma_r = 0.00;
		Sigma_rr = 0.00;
		for (int j = 0; j < p*p; j++)
		{
			kkk = (floor(sqrt(j)) + 0.00);
			Sigma_rr = Sigma_rr + C[j] * Q[j] / (DeltaR + 0.00);
			Sigma_r = Sigma_r + C[j] * Q[j] * (1.0/ (DeltaR + 0.00) +kkk/(Rb+0.0) );
		}
		Potential = Potential + ((Rb+0.0)*(Rb+0.0)) *(Sigma_rr/Distance(x,y,z,DD.x[i],DD.y[i],DD.z[i])-Sigma_r/Distance(x,y,z,D.x[i],D.y[i],D.z[i]))/ (Np+0.00);
		delete Q;
	}
	return -Potential;
}






//This code is used for calculate far field potential by multipole expansion only. It is a method without project.

double CalculateFarFieldPotentialOld(double *C,double x, double y, double z)
{
	double Potential = 0.00;
	double Q[p*p];
	CalculateMultipleExpansion(Q,x,y,z);
	//for (int i = 0; i < p*p; i++)
		//cout << Q[i] << endl;
	for (int i = 0; i < p*p; i++)
	{
		//Potential = Potential + C[i] * P[i](x, y, z);
		Potential = Potential + C[i] * Q[i];
	}
	return Potential;
}
