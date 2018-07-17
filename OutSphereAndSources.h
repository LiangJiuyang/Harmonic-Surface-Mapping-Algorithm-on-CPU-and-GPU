#include"MultipoleExpansion.h"
double Rb = 1.5*sqrt(3.0);
 
void FindAllFieldPointSum(double Radius,double l1,double l2,double l3,Particle *E,Particle *PointBoxSum1)//Radius is the value of Rb. l1 is L1, l2 is L2,l3 is L3. E is set of field charges in center box. N is the number of field charge.
{
        const int N=CenterBoxSourceNumber;

	int x1 = int((Radius - l1 / 2) / l1) + 1;//px=l1,2*l1,...,x1*l1
	int x2 = int((Radius - l2 / 2) / l2) + 1;//py=l2,2*l2,...,x2*l2
	int x3 = int((Radius - l3 / 2) / l3) + 1;//pz=l3,2*l3,...,x3*l3
	int BoxSum = (2 * x1 + 1)*(2 * x2 + 1)*(2 * x3 + 1);//number of box which can wrap sphere Sb
	Particle *PointBoxSum=new Particle[BoxSum*N];
	int l = 0;
	for (int i = -x1; i < x1 + 0.2; i++)
	{
		for (int j = -x2; j < x2 + 0.2; j++)
		{
			for (int k = -x3; k < x3 + 0.2; k++)
			{
				for (int u = 0; u < N; u++)
				{
					PointBoxSum[l + u].q = E[u].q;
					PointBoxSum[l + u].x = E[u].x + i*(l1+0.0);
					PointBoxSum[l + u].y = E[u].y + j*(l2+0.0);
					PointBoxSum[l + u].z = E[u].z + k*(l3+0.0);
				}
				l = l + N;
			}
		}
	}
	//cout << "源电荷总数为=" << BoxSum*N << endl;
	for (int i = 0; i < BoxSum*N; i++)
		PointBoxSum1[i] = PointBoxSum[i];
}

int AdjustAllFieldPointSum(Particle *E,int w,Particle *W)
{
	int q = 0;
	Particle *Y = new Particle[w];
	for (int i = 0; i < w; i++)
	{
		if ((E[i].x*E[i].x+ E[i].y*E[i].y+ E[i].z*E[i].z)<=(Rb*Rb))
		{
			Y[q] = E[i];//Find field charge in Sb
			q++;
		}
	}
	//cout << q << endl;
	for (int i = 0; i < q; i++)
	{
		W[i] = Y[i];
	}
	return q;
}