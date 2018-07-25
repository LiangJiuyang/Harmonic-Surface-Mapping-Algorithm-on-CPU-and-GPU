#include<stdio.h>
#include"CalculateFarField.h"
#include <pthread.h>//线程操作所需头文件
#include <omp.h> 
#include "Eigen/Dense" 
using namespace std;
using namespace Eigen;
struct complex22 { double real, imag; };  //define a complex structure to make the  C code consistent with fortran codes.  Complex*16
extern "C" void lfmm3dparttarg_(int *, int *, int *, double[][3], int *, complex22[], int *, complex22[], double[][3], int *, complex22[], int *, complex22[][3], int *, double[][3], int *, complex22[], int *, complex22[][3]);//State FMM function
//For p>17,we need to set Stack Size to 2M, or will cause stack overflow.
extern "C" void lfmm3dpartself_(int *, int *, int *, double[][3], int *, complex22[], int *, complex22[], double[][3], int *, complex22[], int *, complex22[][3]);//State FMM function

int main()
{
  double startopenmp;  
  double endopenmp;  
  startopenmp = omp_get_wtime(); 

  const int N=ChargeNumber;
  p=25;
  QuizNumber=2*p*p;
  
  //Set Parameter for Fibonacci Integral
  Fp = 610.0;
	F = 987.0;
	Np = F * 2.0 + 2.0;
	deltaz = 2.0 / F; 
  

                  
		clock_t startRead, finishRead;//clock
		startRead = clock();
		double duration;


		Particle FieldCharge[N];
		SetFieldCharge(FieldCharge);//set charge of field source
		QuizPointSum PointSum;// quiz points set (has been constructed in the constructor)

		finishRead = clock();
		printf("Time to read points is:");
		duration = (double)(finishRead - startRead) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
		cout << endl << endl << endl;

            
		cout << "Algorithm include FMM : " << endl;
		cout << "N = " << N << endl;
		cout << "p= " << p << endl;
		cout << "QuizPoints=" << QuizNumber << endl;
		cout << "Np = " << Np << endl;
		cout << "DeltaR=" << DeltaR << endl;

		clock_t startFULL, finishFULL;//clock
		startFULL = clock();

		QuizPointSum QuizSum;//Corresponding quiz points which is moved into center box Omega0 from PointSum that on the circle S0
		QuizSum = FindCorrespondingPoints(PointSum);

		//set L*P fitting matrix A    L=QuizNumber  P=p^2
		 double A[QuizNumber][p*p];//this may cause stack flow

		//double kkk;//use for debug
		double *PointSumMultipleExpansionMatrix = new double[p*p];
		double *QuizSumMultipleExpansionMatrix = new double[p*p];
    double startopenmp33;  
    double endopenmp33;  
    startopenmp33 = omp_get_wtime(); 
		for (int i = 0; i < QuizNumber; i++)
		{
			CalculateMultipleExpansion(PointSumMultipleExpansionMatrix, PointSum.x[i], PointSum.y[i], PointSum.z[i]);
			CalculateMultipleExpansion(QuizSumMultipleExpansionMatrix, QuizSum.x[i], QuizSum.y[i], QuizSum.z[i]);
			for (int j = 0; j < p*p; j++)
			{
				A[i][j] = PointSumMultipleExpansionMatrix[j] - QuizSumMultipleExpansionMatrix[j];
			}
		}
    endopenmp33 = omp_get_wtime();  
    printf("Openmp 统计解方程时间 %f s.\n", endopenmp33-startopenmp33);  
		delete PointSumMultipleExpansionMatrix;
		delete QuizSumMultipleExpansionMatrix;

		//Find all of the period field charges
		int xx1 = int((Rb - L1 / 2) / (L1+0.00)) + 1;//px=L1,2*L1,...,xx1*L1
		int xx2 = int((Rb - L2 / 2) / (L2+0.00)) + 1;//py=L2,2*L2,...,xx2*L2
		int xx3 = int((Rb - L3 / 2) / (L3+0.00)) + 1;//pz=L3,2*L3,...,xx3*L3
		int BBoxSum = (2 * xx1 + 1)*(2 * xx2 + 1)*(2 * xx3 + 1);//number of box which can wrap sphere Sb
    Particle *PP = new Particle[BBoxSum*N];    
		FindAllFieldPointSum(Rb, L1, L2, L3, FieldCharge, PP);//Find all of the field charges. set PP

		int qq = 0;
		Particle *PPP = new Particle[BBoxSum*N];
		qq = AdjustAllFieldPointSum(PP, BBoxSum*N, PPP);//Find field charge in Sb
		Particle *YY = new Particle[qq];
		for (int i = 0; i < qq; i++)
		{
			YY[i] = PPP[i];//Field charges in Sb are set in YY
		}

		clock_t start, finish;//clock
		start = clock();
          
    double startopenmp5;  
    double endopenmp5;  
    startopenmp5 = omp_get_wtime(); 
		/*                                      Begin FMM                                        */
		int iprec = 3;
		int nsource = qq;
		int ifcharge = 1;
		int ifdipole = 0;
		int ifpot = 0;//1
		int iffld = 0;
		int ntarget = QuizNumber + QuizNumber;
		int ifpottarg = 1;
		int iffldtarg = 0;
		//output
		int ier;
                  	
  	//Begin preparing for FMM(Every variable's definition is showed in the end of this page)
		//input(Calculate near field potential)
		//iprec = 1;
		nsource = qq;
		double(*source)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
		complex22 *charge = new complex22[nsource];
		for (int i = 0; i < nsource; i++)
		{
			source[i][0] = YY[i].x;
			source[i][1] = YY[i].y;
			source[i][2] = YY[i].z; 
			charge[i].real = YY[i].q;
			charge[i].imag = 0.0;
		}
		 ifcharge = 1;
		 ifdipole = 0;
		complex22 *dipstr = new complex22[nsource];//do not need to initialize
		double(*dipvec)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));//do not need to initialize
		 ifpot = 0;//1
		 iffld = 0;
		 ntarget = QuizNumber + QuizNumber;
		double(*target)[3] = (double(*)[3])malloc(ntarget * sizeof(double[3]));
		for (int i = 0; i < QuizNumber; i++)
		{
			target[i][0] = QuizSum.x[i];
			target[i][1] = QuizSum.y[i];
			target[i][2] = QuizSum.z[i];

			target[i + QuizNumber][0] = PointSum.x[i];
			target[i + QuizNumber][1] = PointSum.y[i];
			target[i + QuizNumber][2] = PointSum.z[i];
		}
		
		 ifpottarg = 1;
		 iffldtarg = 0;
		//output
		complex22 *pot = new complex22[nsource];
		complex22(*fld)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
		complex22 *pottarg = new complex22[ntarget];
		complex22(*fldtarg)[3] = (complex22(*)[3])malloc(ntarget * sizeof(complex22[3]));
		//End preparing for FMM

		//Begin call FMM procedure in Fortran
		lfmm3dparttarg_(&ier, &iprec, &nsource, source, &ifcharge, charge, &ifdipole, dipstr, dipvec, &ifpot, pot, &iffld, fld, &ntarget, target, &ifpottarg, pottarg, &iffldtarg, fldtarg);//Calculate near field potential for QuizSum
		
		double f[QuizNumber];

		for (int i = 0; i < QuizNumber; i++)
		{
			f[i] = pottarg[i].real - pottarg[i + QuizNumber].real;
		}
                  
    endopenmp5 = omp_get_wtime();  
    printf("Openmp Calculate Matrix A %f s.\n", endopenmp5-startopenmp5);
		finish = clock();
		printf("Time for FMM to calculate matrix A is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);

		double **B;
    B=new double*[QuizNumber];
    for(int i=0;i<QuizNumber;i++)
    B[i]=new double[p*p-1];
    for (int i = 0; i < QuizNumber; i++)
		{
			for (int j = 0; j < p*p - 1; j++)
			{
				if (j < p*p - 1)
				{
					B[i][j] = A[i][j + 1];
				}
			}
		}
        
    /*                        End calculate for multi-pole expansion                            */
    start = clock();
                   
		double *C=new double[p*p];//Storage Solution  
  
	  /*                      Start        Least      Square                         */
     //using Eigen
     double startopenmp22;  
     double endopenmp22;  
     startopenmp22 = omp_get_wtime();
     // Use Eigen instead of mkl               
     MatrixXd MatrixA(QuizNumber,p*p-1);
     VectorXd vB(QuizNumber);
    for (int i = 0; i < QuizNumber; i++)
		{
			for (int j = 0; j < p*p-1; j++)
			{
				if (j < p*p-1)
				{
					MatrixA(i,j) = B[i][j];
				}
			}
                           vB(i)=f[i];
		} 
    MatrixXd MatrixATA=MatrixA.transpose()*MatrixA;
    VectorXd vATB=MatrixA.transpose()*vB;   
    VectorXd vFinal=MatrixATA.colPivHouseholderQr().solve(vATB);    
    for(int i=0;i<p*p-1;i++)
       C[i+1]=vFinal(i);
    C[0]=0;
    endopenmp22 = omp_get_wtime();  
     printf("Openmp 统计解方程时间 %f s.\n", endopenmp22-startopenmp22);  

		/*                      End        Least      Square                         */

		finish = clock();
		printf("Time to solve the LS Least Squares is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
		cout << endl;
                  
		start = clock();

		ProjectPointSum D;//D is points on R'
		ProjectPointSum DD;//DD is corresponding points on R''
                  ProjectPointSum DDD;//DD is corresponding points on R'''
		for (int i = 0; i < Np; i++)
		{
			DD.x[i] = D.x[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);//Rb+0.00 is needed. Please do not use Rb only
			DD.y[i] = D.y[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);
			DD.z[i] = D.z[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);
                          DDD.x[i] = D.x[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);//Rb+0.00 is needed. Please do not use Rb only
			DDD.y[i] = D.y[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);
			DDD.z[i] = D.z[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);
		}


		//Begin FMM(for calculate near field potential between center box and other box)

			//begin FMM output parameters set
		//iprec = 1;
		ntarget = N;
		double(*targetField)[3] = (double(*)[3])malloc(ntarget * sizeof(double[3]));
                  
		for (int i = 0; i < N; i++)
		{
			targetField[i][0] = FieldCharge[i].x;
			targetField[i][1] = FieldCharge[i].y;
			targetField[i][2] = FieldCharge[i].z;
		}
		complex22 *pottargField = new complex22[ntarget];
		complex22(*fldtargField)[3] = (complex22(*)[3])malloc(ntarget * sizeof(complex22[3]));
		//end FMM output parameters set
		//Begin remove the same element in YY(delete N sources in the center box then add their interactions. Or will get number/0=inf) We need to modify source, nsource and charge.
		nsource = nsource - N;
		int Indicate;
		for (int i = 0; i < qq; i++)
		{
			if ((YY[i].x == FieldCharge[0].x) && (YY[i].y == FieldCharge[0].y) && (YY[i].z == FieldCharge[0].z))
				Indicate = i;//YY[Indicate] to YY[Indicete+N-1] is particle in the center box
		}
                                cout<<qq<<"   "<<N<<endl;
		Particle *YYNoCenter = new Particle[qq - N];
		for (int i = 0; i < qq; i++)
		{
			if (i < Indicate)
			{
				YYNoCenter[i] = YY[i];
			}
			else if (i > Indicate + N - 1)
			{
				YYNoCenter[i - N] = YY[i];
			}
		}
		//End remove
 
		//Begin FMM input parameter's set
		nsource = qq - N;
		double(*sourceField)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
		complex22 *chargeField = new complex22[nsource];
		for (int i = 0; i < nsource; i++)
		{
			sourceField[i][0] = YYNoCenter[i].x;
			sourceField[i][1] = YYNoCenter[i].y;
			sourceField[i][2] = YYNoCenter[i].z;
			chargeField[i].real = YYNoCenter[i].q;
			chargeField[i].imag = 0.0;
		}
		complex22 *dipstrField = new complex22[nsource];//do not need to initialize
		double(*dipvecField)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));//do not need to initialize
		complex22 *potField = new complex22[nsource];
		complex22(*fldField)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
		//End FMM input parameter's set

		lfmm3dparttarg_(&ier, &iprec, &nsource, sourceField, &ifcharge, chargeField, &ifdipole, dipstrField, dipvecField, &ifpot, potField, &iffld, fldField, &ntarget, targetField, &ifpottarg, pottargField, &iffldtarg, fldtargField);//Calculate near field potential for PointSum

	//End FMM(for calculate near field potential between center box and other box)
	//Begin FMM( for calculate near field potental in center box)
		//Set Parameter
		//iprec = 1;
		nsource = N;
		double(*sourceCenter)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
		complex22 *chargeCenter = new complex22[nsource];
		for (int i = 0; i < nsource; i++)
		{
			sourceCenter[i][0] = FieldCharge[i].x;
			sourceCenter[i][1] = FieldCharge[i].y;
			sourceCenter[i][2] = FieldCharge[i].z;
			chargeCenter[i].real = FieldCharge[i].q;
			chargeCenter[i].imag = 0.0;
		}
		complex22 *dipstrCenter = new complex22[nsource];//do not need to initialize
		double(*dipvecCenter)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));//do not need to initialize
		ifcharge = 1;
		ifdipole = 0;
		ifpot = 1;//0
		iffld = 0;
		complex22 *potCenter = new complex22[nsource];
		complex22(*fldCenter)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
		//End Set Parameter
		lfmm3dpartself_(&ier, &iprec, &nsource, sourceCenter, &ifcharge, chargeCenter, &ifdipole, dipstrCenter, dipvecCenter, &ifpot, potCenter, &iffld, fldCenter);
		//End FMM(for calculate near field potential in center box)
    finish = clock();
		printf("Time to calculate the near field energy is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
		cout << endl;

                 
     double startopenmp44;  
     double endopenmp44;  
     startopenmp44 = omp_get_wtime(); 

		//Prepare for calculate far field potential
		double *Sigma_r = new double[Np];//Sigma_r is for sigma(r')
		double *Sigma_rr = new double[Np];//Sigma_rr is for sigma(r'')
		double kkk;
		for (int i = 0; i < Np; i++) 
		{
			double *Q = new double[p*p]; 
			CalculateMultipleExpansion(Q, D.x[i], D.y[i], D.z[i]);
			Sigma_r[i] = 0.00;
			Sigma_rr[i] = 0.00;
			for (int j = 0; j < p*p; j++)
			{
				kkk = (floor(sqrt(j)) + 0.00);
				Sigma_rr[i] = Sigma_rr[i] + C[j] * Q[j] / (DeltaR + 0.00);
				Sigma_r[i] = Sigma_r[i] + C[j] * Q[j] * (kkk / (Rb + 0.0));
			}
		}
    endopenmp44 = omp_get_wtime();  
    printf("Openmp 统计解方程时间 %f s.\n", endopenmp44-startopenmp44);  
                
    start=clock();
		//Begin FMM(for calculate far field potential)
		//input(Calculate far field potential)
		nsource = Np;
		double(*sourceFarD)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
		double(*sourceFarDD)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
    double(*sourceFarDDD)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));
		complex22 *chargeFarD = new complex22[nsource];
		complex22 *chargeFarDD = new complex22[nsource];
    complex22 *chargeFarDDD = new complex22[nsource];
		for (int i = 0; i < nsource; i++)
		{
			sourceFarD[i][0] = D.x[i];
			sourceFarD[i][1] = D.y[i];
			sourceFarD[i][2] = D.z[i];
			sourceFarDD[i][0] = DD.x[i];
			sourceFarDD[i][1] = DD.y[i];
			sourceFarDD[i][2] = DD.z[i];
      sourceFarDDD[i][0] = DDD.x[i];
			sourceFarDDD[i][1] = DDD.y[i];
			sourceFarDDD[i][2] = DDD.z[i];
			chargeFarD[i].real = ((Rb + 0.0)*(Rb + 0.0))*D.q[i]*Sigma_r[i] / (4*pi);
			chargeFarD[i].imag = 0.0;
			chargeFarDD[i].real = ((Rb + 0.0)*(Rb + 0.0))*DD.q[i]*Sigma_rr[i] / (4*pi);
			chargeFarDD[i].imag = 0.0;
      chargeFarDDD[i].real = ((Rb + 0.0)*(Rb + 0.0))*DDD.q[i]*Sigma_rr[i] / (4*pi);
			chargeFarDDD[i].imag = 0.0;

		}
		ifcharge = 1;
		ifdipole = 0;
		complex22 *dipstrFar = new complex22[nsource];//do not need to initialize
		double(*dipvecFar)[3] = (double(*)[3])malloc(nsource * sizeof(double[3]));//do not need to initialize
		ifpot = 0;//1
		iffld = 0;
		ntarget = N;
		double(*targetFar)[3] = (double(*)[3])malloc(ntarget * sizeof(double[3]));
		for (int i = 0; i < ntarget; i++)
		{
			targetFar[i][0] = FieldCharge[i].x;
			targetFar[i][1] = FieldCharge[i].y;
			targetFar[i][2] = FieldCharge[i].z;
		}
		ifpottarg = 1;
		iffldtarg = 0;
		//output
		complex22 *potFarD = new complex22[nsource];
		complex22(*fldFarD)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
		complex22 *potFarDD = new complex22[nsource];
		complex22(*fldFarDD)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
    complex22 *potFarDDD = new complex22[nsource];
		complex22(*fldFarDDD)[3] = (complex22(*)[3])malloc(nsource * sizeof(complex22[3]));
		complex22 *pottargFarD = new complex22[ntarget];
		complex22(*fldtargFarD)[3] = (complex22(*)[3])malloc(ntarget * sizeof(complex22[3]));
		complex22 *pottargFarDD = new complex22[ntarget];
		complex22(*fldtargFarDD)[3] = (complex22(*)[3])malloc(ntarget * sizeof(complex22[3]));
    complex22 *pottargFarDDD = new complex22[ntarget];
		complex22(*fldtargFarDDD)[3] = (complex22(*)[3])malloc(ntarget * sizeof(complex22[3]));

		//calling FMM
		lfmm3dparttarg_(&ier, &iprec, &nsource, sourceFarD, &ifcharge, chargeFarD, &ifdipole, dipstrFar, dipvecFar, &ifpot, potFarD, &iffld, fldFarD, &ntarget, targetFar, &ifpottarg, pottargFarD, &iffldtarg, fldtargFarD);//Calculate Far field potential for D
		lfmm3dparttarg_(&ier, &iprec, &nsource, sourceFarDD, &ifcharge, chargeFarDD, &ifdipole, dipstrFar, dipvecFar, &ifpot, potFarDD, &iffld, fldFarDD, &ntarget, targetFar, &ifpottarg, pottargFarDD, &iffldtarg, fldtargFarDD);//Calculate Far field potential for DD
    lfmm3dparttarg_(&ier, &iprec, &nsource, sourceFarDDD, &ifcharge, chargeFarDDD, &ifdipole, dipstrFar, dipvecFar, &ifpot, potFarDDD, &iffld, fldFarDDD, &ntarget, targetFar, &ifpottarg, pottargFarDDD, &iffldtarg, fldtargFarDDD);//Calculate Far field potential for DDD
	//End FMM(for calculate far field potential)

                  
		double tt = 0.00;
    double *Potential=new double[N];
		for (int i = 0; i < N; i++)
		{
		  Potential[i] = pottargField[i].real + potCenter[i].real + pottargFarD[i].real - pottargFarDD[i].real+pottargFarDDD[i].real;
      tt = tt + FieldCharge[i].q * (pottargField[i].real + potCenter[i].real + pottargFarD[i].real - pottargFarDD[i].real+pottargFarDDD[i].real);
		  cout<<Potential[i]<<"  Relative Error is  "<<abs(1-abs(Potential[i])/1.74756459463318219)<<endl;
    }

		finish = clock();
		printf("Time to calculate the final energy is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
		cout << endl;
 
    cout << "The Relative Error is " << abs(1 - ((1.0 + 0.0)*(-1.74756459463318219)) / (tt / (8.0))) << endl;
		cout << "Calculated by simulate is " << setprecision(16) << tt / 8.0 << endl;

		finishFULL = clock();
		printf("Time for all of this prodecure is:");
		duration = (double)(finishFULL - startFULL) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
      
    endopenmp = omp_get_wtime();  
    printf("Openmp 统计总时间 %f s.\n", endopenmp-startopenmp);  
	  return 1;
} 

/*

Input Parameters:

iprec integer :
precision flag. Allowed values are
iprec = −2 for least squares errors < 0.5 100
iprec = −1 for least squares errors < 0.5 10−1
iprec = 0 for least squares errors < 0.5 10−2
iprec = 1 for least squares errors < 0.5 10−3
iprec = 2 for least squares errors < 0.5 10−6
iprec = 3 for least squares errors < 0.5 10−9
iprec = 4 for least squares errors < 0.5 10−12
iprec = 5 for least squares errors < 0.5 10−15

nsource integer :
number of sources

source(3,nsources) real *8 :
sources(k,j) is the kth component of the jth source in R3

ifcharge integer :
charge flag. If icharge = 1, then include the effect of the charge sources. Otherwise,
omit.

charge(nsources) complex *16 :
charge(j) is the strength of the jth charge (qj in the formula (1)).

ifdipole integer :
dipole flag. If idipole = 1, then include the effect of the dipole sources. Otherwise,
omit.

dipstr(nsources) complex *16 :

dipstr(j) is the strength of the jth dipole (pj in the formula (1)).

dipvec(3,nsources) real *8 :

dipvec(k,j) is the kth component of the orientation vector of the jth dipole (nj in the
formula (1)).

ifpot integer :
potential flag. If ifpot = 1, the potential is computed. Otherwise, it is not.

iffld integer :
field (gradient) flag. If iffld = 1 the gradient of the potential is computed. Other wise, it is not.

ntarget integer :
number of targets

target(3,ntarget) real *8 :

target(k,j) is the kth component of the jth target in R3

ifpottarg integer :
target potential flag. If ifpottarg = 1, the potential is computed. Otherwise, it is
not.

iffldtarg integer :
target field (gradient) flag. If iffldtarg = 1 the gradient of the potential is computed.
Otherwise, it is not.
Unused arrays do not need to be allocated in full. Thus, if ifcharge = 0, charge
can be dimensioned as a (complex) scalar. If ifdipole = 0, dipstr can be di mensioned as a complex scalar and dipvec can be dimensioned in the calling
program as dipvec(3) - BUT NOT dipvec(1).

Output Parameters:

ier integer :
Error return codes.
ier = 0: Successful completion of code.
ier = 4: failure to allocate memory for oct-tree
ier = 8: failure to allocate memory for FMM workspaces
ier = 16: failure to allocate meory for multipole/local expansions

pot(nsources) complex *16 :
pot(i) is the potential at the ith source

fld(3,nsources) complex *16 :
fld(k,i) is the kth component of the field (-gradient of the potential) at the ith source

pottarg(ntarget) complex *16 :
pottarg(i) is the potential at the ith target

fldtarg(3,ntarget) complex *16 :
fldtarg(k,i) is the kth component of the field (-gradient of the potential) at the ith
target

Note that the charge, dipstr, pot, fld, pottarg, fldtarg, arrays must be declared
and passed as complex arrays (even if the charge and dipole strengths are real)

*/
