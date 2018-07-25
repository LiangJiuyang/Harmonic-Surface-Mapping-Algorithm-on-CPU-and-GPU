#include<stdio.h>
#include"CalculateFarField.h"
#include<cuda_runtime.h>
#include "cublas_v2.h"
#include <pthread.h>//head file needed for thread operation
using namespace std;

//MKL macro definition
#ifndef lapack_int
_#define lapack_int MKL_INT
#endif
#ifndef lapack_logical
_#define lapack_logical lapack_int
#endif

extern void line_fun(double **a,double *b,double *ppp);
int BlockNum;//Number of blocks
int ThreadNum;//Number of threads

/*
    Set number of GPU. Note that this code is based on two GPU whose cudaStatus number are 0 and 2. 
    User need to find your own GPU's cudaStatus number.
    For example, one have a GPU which cudaStatus 0, all of the code "cudaStatus = cudaSetDevice(2)" need to be changed to "cudaStatus = cudaSetDevice(0)".  
*/
cudaError_t cudaStatus;
 
//bodyBody interaction qj/|ri-rj|
__device__ double bodyBodyInteraction(double4 bi, double4 bj, double ai)
{
	double3 r;

	//r_ij 
	r.x = bj.x - bi.x;
	r.y = bj.y - bi.y;
	r.z = bj.z - bi.z;

	double distSqr = r.x*r.x + r.y*r.y + r.z*r.z;

	double invDistCube = (bj.w) / sqrt(distSqr);

	if (r.x == 0 && r.y == 0 && r.z == 0)
	  {
		  invDistCube = 0.0;
	  }
	ai = ai + invDistCube;
	return ai;
}

//tile calculation(sum of body-body interactions in a tile)
__device__ double tile_calculation(double4 myPosition, double accel)
{
	int i;
	extern __shared__ double4 shPosition[];
	for (i = 0; i<blockDim.x; i++)
	{
		accel = bodyBodyInteraction(myPosition, shPosition[i], accel);
	}
	return accel;
}

//sum of body-body interactions in all of the tiles
__global__ void calculate_forces(double4 *devX, double4 *devY, double *devA, int *Number)
{
	extern __shared__ double4 shPosition[];
	double4 myPosition;
	int i, tile;
	double acc = 0.0;
	int gtid = blockIdx.x*blockDim.x + threadIdx.x;

	if (gtid<Number[0])
	{
		myPosition = devX[gtid];
	}

	for (i = 0, tile = 0; i<Number[1]; i += blockDim.x, tile++)
	{
		int idx = tile*blockDim.x + threadIdx.x;

		if (idx<Number[1])
		{
			shPosition[threadIdx.x] = devY[idx];
		}
		else if (idx >= Number[1])
		{
			shPosition[threadIdx.x].x = 0.00;
			shPosition[threadIdx.x].y = 0.00;
			shPosition[threadIdx.x].z = 0.00;
			shPosition[threadIdx.x].w = 0.00;
		}

		__syncthreads();

		if (gtid<Number[0])
		{
			acc = tile_calculation(myPosition, acc);
		}

		__syncthreads();
	}

	devA[gtid] = acc;

}

//Energy calculation(these three components are many same as body-body interaction) qiqj/|ri-rj|
__device__ double bodyBodyInteraction1(double4 bi, double4 bj, double ai)
{
	double3 r;

	//r_ij 
	r.x = bj.x - bi.x;
	r.y = bj.y - bi.y;
	r.z = bj.z - bi.z;

	double distSqr = r.x*r.x + r.y*r.y + r.z*r.z;
	double invDistCube = (bi.w+0.00)*(bj.w + 0.00) / sqrt(distSqr);
	if (r.x == 0 && r.y == 0 && r.z == 0)
	  {
		  invDistCube = 0.0;
	  }
	ai = ai + invDistCube;
	return ai;
}

__device__ double tile_calculation1(double4 myPosition, double accel)
{
	int i;
	extern __shared__ double4 shPosition[];
	for (i = 0; i<blockDim.x; i++)
	  {
	   	accel = bodyBodyInteraction1(myPosition, shPosition[i], accel);
  	}
	return accel;
}

__global__ void calculate_forces1(double4 *devX, double4 *devY, double *devA, int *Number)
{
	extern __shared__ double4 shPosition[];
	double4 myPosition;
	int i, tile;
	double acc = 0.0;
	int gtid = blockIdx.x*blockDim.x + threadIdx.x;
	if (gtid<Number[0])
	{
		myPosition = devX[gtid];
	}
	else if(gtid>=Number[0])
	{
	        myPosition.x = 1.00;
			myPosition.y = 1.00;
			myPosition.z = 1.00;
			myPosition.w = 0.00;
	}

	for (i = 0, tile = 0; i<Number[1]; i += blockDim.x, tile++)
	{
		int idx = tile*blockDim.x + threadIdx.x;

		if (idx<Number[1])
		{
			shPosition[threadIdx.x] = devY[idx];
		}
		else if (idx >= Number[1])
		{
			shPosition[threadIdx.x].x = 0.00;
			shPosition[threadIdx.x].y = 0.00;
			shPosition[threadIdx.x].z = 0.00;
			shPosition[threadIdx.x].w = 0.00;
		}

		__syncthreads();
	}

	devA[gtid] = acc;	
}

//Calculate factorial on GPU 
__device__ double facGpu(double ttt)
{
    double s;
    if (abs(ttt - 1)<0.001 || abs(ttt)<0.001)
		s = 1.0;
	else
	{ 
		s = ttt*facGpu(ttt - 1) + 0.00;
	}
	return s;
}

//Calculate multipole expansion on GPU
__global__ void Calculate_MultipoleExpansion(double *ExpansionMatrix,double4 *PointPosition,int *QuizPointNumber,int *ExpansionNumber)
{
           double4 myPosition;
           int gtid = blockIdx.x*blockDim.x + threadIdx.x;
           myPosition=PointPosition[gtid];
           double x,y,z;
           x=myPosition.x;
           y=myPosition.y;
           z=myPosition.z;
           if(gtid<QuizPointNumber[0]) 
       {             
	ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+0] = 1.0;
	ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+1] = y / 2;
	ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+2] =(-z);
	ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+3] =(-x / 2);
	int t = 4;
	int m, n,i;
	for (i = 2; i < ExpansionNumber[0]; i++) 
	{
		while (t<(i+1)*(i+1))
		{
			m = t - i - i*i;
			n = i; 
			if (m==-n)
			{
				ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] =(y*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n - 1] - x*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n - 2 * n + 1]) / (2 * n + 0.00);
			}
			else if (m==n)
			{
				ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] = (-(x*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n - 1] + y*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n - 2 * n + 1]) / (2 * n + 0.00));
			}
			else if (n - abs(m) == 1)
			{
				ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] = (-z)*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n-n+m];
			}      
			else if ((n-abs(m))>1)    
			{
				ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] = (-((2 * n - 1.0)*z*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n-n+m] + (x*x + y*y + z*z)*ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+n*n-3*n+m+2]) / ((n - abs(m)+0.0)*(n + abs(m)+0.0)));
			}
			t++;
		}
	}

	t = 0;
	for (i = 0; i < ExpansionNumber[0]; i++)
	{
		while (t < (i + 1)*(i + 1))
		{
			m = t - i - i*i;
			n = i;
			ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] = ExpansionMatrix[gtid*ExpansionNumber[0]*ExpansionNumber[0]+t] * sqrt((facGpu(n - m+0.0))*facGpu(n + m+0.0));
			t++;
		}
	}
       }
}

//Output device information
void printDeviceProp(const cudaDeviceProp &prop)
{
    printf("Device Name : %s.\n", prop.name);
    printf("totalGlobalMem : %d.\n", prop.totalGlobalMem);
    printf("sharedMemPerBlock : %d.\n", prop.sharedMemPerBlock);
    printf("regsPerBlock : %d.\n", prop.regsPerBlock);
    printf("warpSize : %d.\n", prop.warpSize);
    printf("memPitch : %d.\n", prop.memPitch);
    printf("maxThreadsPerBlock : %d.\n", prop.maxThreadsPerBlock);
    printf("maxThreadsDim[0 - 2] : %d %d %d.\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
    printf("maxGridSize[0 - 2] : %d %d %d.\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    printf("totalConstMem : %d.\n", prop.totalConstMem);
    printf("major.minor : %d.%d.\n", prop.major, prop.minor);
    printf("clockRate : %d.\n", prop.clockRate);
    printf("textureAlignment : %d.\n", prop.textureAlignment);
    printf("deviceOverlap : %d.\n", prop.deviceOverlap);
    printf("multiProcessorCount : %d.\n", prop.multiProcessorCount);
}

bool InitCUDA()
{
    //used to count the device numbers
    int count;

    // get the cuda device count
    cudaGetDeviceCount(&count);
    if (count == 0) {
        fprintf(stderr, "There is no device.\n");
        return false;
    }

    // find the device >= 1.X

    int i;
    for (i = 0; i < count; ++i) {
        cudaDeviceProp prop;
        if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
            if (prop.major >= 0) {
                printDeviceProp(prop);
                break;
            }
        }
    }

    // if can't find the device
    if (i == count) {
        fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
        return false;
    }

    // set cuda device
    cudaSetDevice(i);
    return true;
}



int main()
{
  //Initial GPU
  if (!InitCUDA()) 
     {
       return 0;
     }
  cudaStatus = cudaSetDevice(2);
  int *kkkkl;
	cudaMalloc((void**)&kkkkl, sizeof(int));
  cudaFree(kkkkl); 
       
	cudaError_t cudaStatus;  
  int num;   
  cudaStatus = cudaGetDeviceCount(&num);
	cout<<"Num = "<<num<<endl; 
	   
  //Set Fibonacci Integral's parameter   
  Fp = 610.0;
  F = 987.0;
	Np = F * 2.0 + 2.0;
	deltaz = 2.0 / F;
	
  //HSMA Parameters' set
  p=25;//number of multiple expansion's truncate
  Rb=1.5*sqrt(3);//Rb is Rs in the related paper     
	QuizNumber=2*p*p;//number of test point on the surface
  const int N=CenterBoxSourceNumber;//number of particles in center box
		
  double duration;
  clock_t start, finish;//clock
                  
  
		Particle FieldCharge[N];
		SetFieldCharge(FieldCharge);//set charge of field source
		QuizPointSum PointSum;// quiz points set (has been constructed in the constructor)


		cout << "Algorithm include FMM : " << endl;
		cout << "N = " << N << endl;
		cout << "p= " << p << endl;
		cout << "Rb = " << (Rb+0.00)/sqrt(3)<<"*sqrt(3)" << endl;
		cout << "QuizPoints=" << QuizNumber << endl;
		cout << "Np = " << Np << endl;
		cout << "DeltaR=" << DeltaR << endl; 
		 
		clock_t startFULL, finishFULL;//clock
		startFULL = clock();
    start=clock();
		
    QuizPointSum QuizSum;//Corresponding quiz points which is moved into center box Omega0 from PointSum that on the circle S0
		QuizSum = FindCorrespondingPoints(PointSum);

		//Find all of the period field charges
		int xx1 = int((Rb - L1 / 2) / L1) + 1;//px=L1,2*L1,...,xx1*L1
		int xx2 = int((Rb - L2 / 2) / L2) + 1;//py=L2,2*L2,...,xx2*L2
		int xx3 = int((Rb - L3 / 2) / L3) + 1;//pz=L3,2*L3,...,xx3*L3
		int BBoxSum = (2 * xx1 + 1)*(2 * xx2 + 1)*(2 * xx3 + 1);//number of box which can wrap sphere Sb
		Particle *PP = new Particle[BBoxSum*N];
		FindAllFieldPointSum(Rb, L1, L2, L3, FieldCharge, PP);//Find all of the field charges. set PP
 
                   
		int qq = 0;
		Particle *PPP = new Particle[BBoxSum*N];
		qq = AdjustAllFieldPointSum(PP, BBoxSum*N, PPP);//Find field charge in Rs
		Particle *YY = new Particle[qq];
		for (int i = 0; i < qq; i++)
		{
			YY[i] = PPP[i];//Field charges in Rs are set in YY
		}

		double *f=new double[QuizNumber];//store right-hand term
  
    finish = clock();
		printf("Time for CPU to initial set is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);

               

/*      begin    cuda       */
    start=clock();
		//Set GPU's Number 
    // Choose which GPU to run on, change this on a multi-GPU system.  
    cudaStatus = cudaSetDevice(1);  
    if (cudaStatus != cudaSuccess) {  
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");    
        } 

		/*   The First GPU   */
		//CPUset
		double(*QuizPointCopy)[4] = (double(*)[4])malloc(QuizNumber * 2 * sizeof(double[4]));//change structure for gpu
		for (int i = 0; i < QuizNumber; i++)
		{
			QuizPointCopy[i][0] = QuizSum.x[i];
			QuizPointCopy[i][1] = QuizSum.y[i];
			QuizPointCopy[i][2] = QuizSum.z[i];
			QuizPointCopy[i][3] = 0.00;
			QuizPointCopy[i + QuizNumber][0] = PointSum.x[i];
			QuizPointCopy[i + QuizNumber][1] = PointSum.y[i];
			QuizPointCopy[i + QuizNumber][2] = PointSum.z[i];
			QuizPointCopy[i + QuizNumber][3] = 0.00;
		}
		//GPUset
		double4 *PositionForQuizPoint;
		double4 *Source,*Source1;//We will use Source twice. So that it do not need to free quickly.
 
		cudaMalloc((void**)&PositionForQuizPoint, sizeof(double4)*QuizNumber * 2);
		//cudaMalloc((void**)&PositionForPoint, sizeof(double4)*QuizNumber);
		cudaMalloc((void**)&Source, sizeof(double4)*qq);
		cudaMalloc((void**)&Source1, sizeof(double4)*qq/2);
		cudaMemcpy(PositionForQuizPoint, QuizPointCopy, sizeof(double4)*QuizNumber * 2, cudaMemcpyHostToDevice);
		cudaMemcpy(Source, YY, sizeof(double4)*qq, cudaMemcpyHostToDevice);
		cudaMemcpy(Source1, YY, sizeof(double4)*qq/2, cudaMemcpyHostToDevice);


		int number[2];
		int *Number;
		number[0] = QuizNumber * 2;
		number[1] = qq/2;
		cudaMalloc((void**)&Number, sizeof(int) * 2);
		cudaMemcpy(Number, number, sizeof(int) * 2, cudaMemcpyHostToDevice);

		double *GPUQuizPointPotential;
		cudaMalloc((void**)&GPUQuizPointPotential, sizeof(double)*QuizNumber * 2);

		BlockNum = QuizNumber*2/128+1;
		ThreadNum = 128;

		//Calculate
		calculate_forces << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (PositionForQuizPoint, Source1, GPUQuizPointPotential, Number);


		//Copy from GPU to CPU
		double GPUQUIZPOINTPOTENTIAL[QuizNumber * 2];

    //Initial the second GPU and substract it from the total time
    clock_t cutstart,cutfinish;
    cutstart=clock();
    cudaStatus = cudaSetDevice(0); 
  	int *kkkl;
	  cudaMalloc((void**)&kkkl, sizeof(int));
	  cudaFree(kkkl);
    cutfinish=clock();
    double cut1;
    cut1=(double)(cutfinish - cutstart) / CLOCKS_PER_SEC;
    cout<<"cut1 = "<<cut1<<endl;

		/*   The Second GPU   */
		double4 *PositionForQuizPoint2;
		double4 *Source2;//We will use Source twice. So that it do not need to free quickly.
		cudaMalloc((void**)&PositionForQuizPoint2, sizeof(double4)*QuizNumber * 2);
		cudaMalloc((void**)&Source2, sizeof(double4)*(qq-qq/2));
		cudaMemcpy(PositionForQuizPoint2, QuizPointCopy, sizeof(double4)*QuizNumber * 2, cudaMemcpyHostToDevice);
		cudaMemcpy(Source2, &YY[qq/2], sizeof(double4)*(qq-qq/2), cudaMemcpyHostToDevice);

                  
		int *Number2;
		number[0] = QuizNumber * 2;
		number[1] = qq-qq/2;
		cudaMalloc((void**)&Number2, sizeof(int) * 2);
		cudaMemcpy(Number2, number, sizeof(int) * 2, cudaMemcpyHostToDevice);
		double *GPUQuizPointPotential2;
		cudaMalloc((void**)&GPUQuizPointPotential2, sizeof(double)*QuizNumber * 2);
		BlockNum = QuizNumber*2/128+1;
		ThreadNum = 128;
		//Calculate
		calculate_forces << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (PositionForQuizPoint2, Source2, GPUQuizPointPotential2, Number2);

    //Copy from GPU to CPU
		double GPUQUIZPOINTPOTENTIAL2[QuizNumber * 2];
		cudaMemcpy(&GPUQUIZPOINTPOTENTIAL, GPUQuizPointPotential, sizeof(double)*QuizNumber * 2, cudaMemcpyDeviceToHost);
		cudaMemcpy(&GPUQUIZPOINTPOTENTIAL2, GPUQuizPointPotential2, sizeof(double)*QuizNumber * 2, cudaMemcpyDeviceToHost);


		
/*                    Prepare for next use                                   */
    double4 *SourceSecond;//We will use Source twice. So that it do not need to free quickly.
    cudaMalloc((void**)&SourceSecond, sizeof(double4)*qq); 
    cudaMemcpy(SourceSecond, YY, sizeof(double4)*qq, cudaMemcpyHostToDevice);
    double4 *PositionFieldSecond;//We will use it twice
		cudaMalloc((void**)&PositionFieldSecond, sizeof(double4)*N/2);
		cudaMemcpy(PositionFieldSecond, &FieldCharge[N/2], sizeof(double4)*N/2, cudaMemcpyHostToDevice);
    int *NumberSecond;
		cudaMalloc((void**)&NumberSecond, sizeof(int)*2);
		number[0] = qq;
		number[1] = N/2;
		cudaMemcpy(NumberSecond, number, sizeof(int) * 2, cudaMemcpyHostToDevice);
/*                    End prepare                                      */


		//Free
		cudaFree(PositionForQuizPoint); cudaFree(GPUQuizPointPotential);
/*       end      cuda       */

  
		for (int i = 0; i < QuizNumber; i++)
		{
			f[i] = GPUQUIZPOINTPOTENTIAL[i] - GPUQUIZPOINTPOTENTIAL[i + QuizNumber]+GPUQUIZPOINTPOTENTIAL2[i] - GPUQUIZPOINTPOTENTIAL2[i + QuizNumber];
		}


    finish = clock();
		printf("Time for GPU to calculate the right part is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration-cut1);


    start=clock();
    cudaStatus = cudaSetDevice(1);
    double *QuizSumExpansionMatrix;
    cudaMalloc((void **)&QuizSumExpansionMatrix,sizeof(double)*QuizNumber*p*p);
    double4 *QuizSumPosition;              
    cudaMalloc((void **)&QuizSumPosition,sizeof(double4)*QuizNumber);
    cudaMemcpy(QuizSumPosition, QuizPointCopy, sizeof(double4) * QuizNumber, cudaMemcpyHostToDevice);
    int  *QuizNumberCuda1;
    int  *ExpansionNumberCuda1; 
    cudaMalloc((void **)&QuizNumberCuda1,sizeof(int));                                
    cudaMalloc((void **)&ExpansionNumberCuda1,sizeof(int));
    int *UM=new int[2];
    UM[0]=QuizNumber;
    UM[1]=p;
    cudaMemcpy(QuizNumberCuda1,&UM[0],sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(ExpansionNumberCuda1,&UM[1],sizeof(int),cudaMemcpyHostToDevice);     
    BlockNum=QuizNumber/128+1;
    ThreadNum=128;
    Calculate_MultipoleExpansion<<<BlockNum,ThreadNum,0>>>(QuizSumExpansionMatrix,QuizSumPosition,QuizNumberCuda1,ExpansionNumberCuda1);

    cudaStatus = cudaSetDevice(0);
    double *PointSumExpansionMatrix;
    cudaMalloc((void **)&PointSumExpansionMatrix,sizeof(double)*QuizNumber*p*p);
    double4 *PointSumPosition;
    cudaMalloc((void **)&PointSumPosition,sizeof(double4)*QuizNumber); 
    cudaMemcpy(PointSumPosition, &QuizPointCopy[QuizNumber], sizeof(double4) * QuizNumber, cudaMemcpyHostToDevice);                                                  
    int *QuizNumberCuda2,*ExpansionNumberCuda2;
    cudaMalloc((void **)&QuizNumberCuda2,sizeof(int)); 
    cudaMalloc((void **)&ExpansionNumberCuda2,sizeof(int));
    cudaMemcpy(QuizNumberCuda2,&UM[0],sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(ExpansionNumberCuda2,&UM[1],sizeof(int),cudaMemcpyHostToDevice);                   
    Calculate_MultipoleExpansion<<<BlockNum,ThreadNum,0>>>(PointSumExpansionMatrix,PointSumPosition,QuizNumberCuda2,ExpansionNumberCuda2);
                  
                    
    //Calculate Projection in cpu
		ProjectPointSum D;//D is points on R'
		ProjectPointSum DD;//DD is corresponding points on R''
    ProjectPointSum DDD;//DD is corresponding points on R'''
    double4 *D4=new double4[Np];
		for (int i = 0; i < Np; i++)
		{
			DD.x[i] = D.x[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);//Rb+0.00 is needed. Please do not use Rb only
			DD.y[i] = D.y[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);
			DD.z[i] = D.z[i] * (Rb + DeltaR/2.0 + 0.0) / (Rb + 0.0);
      DDD.x[i] = D.x[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);//Rb+0.00 is needed. Please do not use Rb only
			DDD.y[i] = D.y[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);
			DDD.z[i] = D.z[i] * (Rb - DeltaR/2.0 + 0.0) / (Rb + 0.0);
      D4[i].x=D.x[i];
      D4[i].y=D.y[i];
      D4[i].z=D.z[i];
      D4[i].w=0.00;
		} 
                  
		//Prepare for calculate far field potential
		double *Sigma_r = new double[Np];//Sigma_r is for sigma(r')
		double *Sigma_rr = new double[Np];//Sigma_rr is for sigma(r'')
		double kkk;
  //End Calculate in cpu    


    double A[QuizNumber][p*p];
    double AA[QuizNumber][p*p];
    cudaStatus = cudaSetDevice(1);
    cudaMemcpy(&A, QuizSumExpansionMatrix, sizeof(double)*QuizNumber*p*p, cudaMemcpyDeviceToHost);                       
    cudaStatus = cudaSetDevice(0);
    cudaMemcpy(&AA, PointSumExpansionMatrix, sizeof(double)*QuizNumber*p*p, cudaMemcpyDeviceToHost);       
    cudaFree(QuizSumExpansionMatrix);
    cudaFree(PointSumExpansionMatrix);
    cudaFree(QuizSumPosition);cudaFree(QuizNumberCuda1);cudaFree(ExpansionNumberCuda1);cudaFree(PointSumPosition);cudaFree(QuizNumberCuda2);cudaFree(ExpansionNumberCuda2);
                                  
    for(int i=0;i<QuizNumber;i++)
    for(int j=0;j<p*p;j++)
      {
        A[i][j]=AA[i][j]-A[i][j];
      }
  
    finish = clock();
		printf("Time for CPU to calculate multi-pole expansion is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);            

    start=clock();   
                   
    cudaStatus = cudaSetDevice(1);
    double *DSumExpansionMatrix;
    cudaMalloc((void **)&DSumExpansionMatrix,sizeof(double)*(Np/2)*p*p);
    double4 *DSumPosition;              
    cudaMalloc((void **)&DSumPosition,sizeof(double4)*(Np/2));
    cudaMemcpy(DSumPosition, D4, sizeof(double4) * (Np/2), cudaMemcpyHostToDevice);
    int  *DNumberCuda1;
    int  *DExpansionNumberCuda1; 
    cudaMalloc((void **)&DNumberCuda1,sizeof(int));                                
    cudaMalloc((void **)&DExpansionNumberCuda1,sizeof(int));
    UM[0]=Np/2;
    UM[1]=p;
    cudaMemcpy(DNumberCuda1,&UM[0],sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(DExpansionNumberCuda1,&UM[1],sizeof(int),cudaMemcpyHostToDevice);     
    BlockNum=(Np/2)/64+1;
    ThreadNum=64;
    Calculate_MultipoleExpansion<<<BlockNum,ThreadNum,0>>>(DSumExpansionMatrix,DSumPosition,DNumberCuda1,DExpansionNumberCuda1);
                     
    cudaStatus = cudaSetDevice(0); 
    BlockNum=(Np-Np/2)/64+1;
    ThreadNum=64;
    double *D2SumExpansionMatrix;
    cudaMalloc((void **)&D2SumExpansionMatrix,sizeof(double)*(Np-Np/2)*p*p);
    double4 *D2SumPosition;   
    cudaMalloc((void **)&D2SumPosition,sizeof(double4)*(Np-Np/2));
    cudaMemcpy(D2SumPosition, &D4[Np/2], sizeof(double4) * (Np-Np/2), cudaMemcpyHostToDevice);                                                    
    int *D2NumberCuda2,*D2ExpansionNumberCuda2;
    cudaMalloc((void **)&D2NumberCuda2,sizeof(int)); 
    cudaMalloc((void **)&D2ExpansionNumberCuda2,sizeof(int));   
    UM[0]=Np-Np/2;
    UM[1]=p;
    cudaMemcpy(D2NumberCuda2,&UM[0],sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(D2ExpansionNumberCuda2,&UM[1],sizeof(int),cudaMemcpyHostToDevice);                   
    Calculate_MultipoleExpansion<<<BlockNum,ThreadNum,0>>>(D2SumExpansionMatrix,D2SumPosition,D2NumberCuda2,D2ExpansionNumberCuda2);

    double AB[Np][p*p];
    cudaStatus = cudaSetDevice(1);
    cudaMemcpy(&AB, DSumExpansionMatrix, sizeof(double)*(Np/2)*p*p, cudaMemcpyDeviceToHost);            
    cudaStatus = cudaSetDevice(0);
    cudaMemcpy(&AB[Np/2][0], D2SumExpansionMatrix, sizeof(double)*(Np-Np/2)*p*p, cudaMemcpyDeviceToHost);        
    cudaFree(DSumExpansionMatrix);   cudaFree(DSumPosition);  cudaFree(DNumberCuda1);  cudaFree(DExpansionNumberCuda1); 
    cudaFree(D2SumExpansionMatrix);   cudaFree(D2SumPosition);  cudaFree(D2NumberCuda2);  cudaFree(D2ExpansionNumberCuda2);                  

    finish = clock();
		printf("Time for CPU to calculate projection is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);
  
    clock_t start1, finish1;//clock
		start1 = clock();	
                  
		//Calculate Near Field Potential
/*                begin         GPU                      */
		cudaStatus = cudaSetDevice(1);  
    /*      The First GPU     */
		//CPUset
		//We have copied YY to GPU.(Source) And we can copy FieldCharge to GPU directly.
		//GPUset
                                 
		double4 *PositionField;//We will use it twice
		cudaMalloc((void**)&PositionField, sizeof(double4)*N/2);
		cudaMemcpy(PositionField, FieldCharge, sizeof(double4)*N/2, cudaMemcpyHostToDevice);
        
		number[0] = qq;
		number[1] = N/2;
		cudaMemcpy(Number, number, sizeof(int) * 2, cudaMemcpyHostToDevice);

		double *GPUFieldPotential;//We will use it twice
		cudaMalloc((void**)&GPUFieldPotential, sizeof(double)*qq);

		BlockNum = qq / 256 + 1;
		ThreadNum = 256;
  
		//Calculate
		calculate_forces1 << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (Source,PositionField,GPUFieldPotential, Number);

		 //Copy from GPU to CPU
		double GPUFIELDPOTENTIAL[qq];//If N=216000, we need to set stack 81920Kbts.
		clock_t time_use1[BlockNum*2];
                  	
		/*      The second gpu      */
    cudaStatus = cudaSetDevice(0); 
		//CPUset
		//We have copied YY to GPU.(Source) And we can copy FieldCharge to GPU directly.
		//GPUset

		double *GPUFieldPotentialSecond;//We will use it twice
		cudaMalloc((void**)&GPUFieldPotentialSecond, sizeof(double)*qq);
                  
		BlockNum = qq / 256 + 1;
		ThreadNum = 256;

		//Calculate
		calculate_forces1 << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (SourceSecond,PositionFieldSecond,GPUFieldPotentialSecond, NumberSecond);
 
		 //Copy from GPU to CPU
		double GPUFIELDPOTENTIALSECOND[qq];//If N=216000, we need to set stack 81920Kbts.
		
/*            end    GPU                */
                   
/*                       Begin calculate for multi-pole expansion                          */
                 
    start=clock();
		//A is the matrix of Quiz - Point
    //Set augmented matrix(used for LS least Squares)
		double(*B)[p*p - 1] = (double(*)[p*p - 1])malloc(QuizNumber * sizeof(double[p*p - 1]));
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

/*                           Solve      Least      Square                             */
    //Using MKL for accerate 
                 
    double *MatrixAT,*MatrixATA;
    int rowAT,columnATA,columnAT; 
    double alpha,beta;
    rowAT=p*p-1;columnAT=QuizNumber;columnATA=p*p-1;
    alpha=1.0;beta=0.00;
    MatrixAT = (double *)mkl_malloc( rowAT*columnAT*sizeof( double ), 64 );
    MatrixATA=(double *)mkl_malloc(rowAT*columnATA*sizeof(double),64);
    if (MatrixAT == NULL || MatrixATA == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(MatrixAT);
      mkl_free(MatrixATA);
      return 1;
    }
    for(int i=0;i<rowAT;i++)
       for(int j=0;j<columnAT;j++)
           {
              MatrixAT[i*columnAT+j]=B[j][i];
            }
    for (int i = 0; i < (rowAT*columnATA); i++) {
         MatrixATA[i] = 0.0;
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowAT, columnATA, columnAT, alpha, MatrixAT, columnAT, MatrixAT, columnAT, beta, MatrixATA, columnATA);
    int InfoHelp,*VectorHelp;
    VectorHelp=(int *)mkl_malloc(columnATA*sizeof(int),64); 
    for(int i=0;i<columnATA;i++)
        VectorHelp[i]=0;
    InfoHelp=LAPACKE_dgetrf(CblasRowMajor,columnATA,columnATA,MatrixATA,columnATA,VectorHelp);
    InfoHelp=LAPACKE_dgetri(CblasRowMajor,columnATA,MatrixATA,columnATA,VectorHelp); 
    double *BB,*ATB,*INV_ATA_ATB;
    BB = (double *)mkl_malloc( columnAT*1*sizeof( double ), 64 );
    ATB = (double *)mkl_malloc( columnATA*1*sizeof( double ), 64 );
    INV_ATA_ATB = (double *)mkl_malloc( columnATA*1*sizeof( double ), 64 );
    for(int i=0;i<columnAT;i++)
        {
           BB[i]=f[i];
        }
    for(int i=0;i<columnATA;i++)
        {
           ATB[i]=0.00;
           INV_ATA_ATB[i]=0.00;
        }
                                                               
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rowAT,columnAT, alpha, MatrixAT, columnAT, BB, 1, beta, ATB, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,columnATA,columnATA, alpha, MatrixATA, columnATA, ATB, 1, beta, INV_ATA_ATB, 1);
                  
		/*                      End        Least      Square                         */
		for (int i = 1; i < p*p; i++)
			C[i] = INV_ATA_ATB[i-1];
      C[0] = 0;
                                                      
		finish = clock();
		printf("Time to solve the LS Least Squares is:");
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);

    //Calculate Parameter for Harmonic Mapping
		for (int i = 0; i < Np; i++)
		{
			Sigma_r[i] = 0.00;
			Sigma_rr[i] = 0.00;
			for (int j = 0; j < p*p; j++)
			{
				kkk = (floor(sqrt(j)) + 0.00);
				Sigma_rr[i] = Sigma_rr[i] + C[j] * AB[i][j] / (DeltaR + 0.00);
				Sigma_r[i] = Sigma_r[i] + C[j] * AB[i][j] * (kkk / (Rb + 0.0));
			}
		 }
                 	
        

		/*       begin CUDA            */
    cudaStatus = cudaSetDevice(1);  
    double4 *PositionField1;//We will use it twice
		cudaMalloc((void**)&PositionField1, sizeof(double4)*N);
		cudaMemcpy(PositionField1, FieldCharge, sizeof(double4)*N, cudaMemcpyHostToDevice);

		int *Number1;
		cudaMalloc((void**)&Number1, sizeof(int) * 2);

		//CPUset
		double(*DCopy)[4] = (double(*)[4])malloc(Np * sizeof(double[4]));
    double(*DDCopy)[4] = (double(*)[4])malloc(Np * sizeof(double[4]));
    double(*DDDCopy)[4] = (double(*)[4])malloc(Np * sizeof(double[4]));
		for (int i = 0; i < Np; i++)
		{ 
			DCopy[i][0] = D.x[i];
			DCopy[i][1] = D.y[i];
			DCopy[i][2] = D.z[i];
			DCopy[i][3] = ((Rb + 0.0)*(Rb + 0.0))*D.q[i] * Sigma_r[i] / (4 * pi);
			DDCopy[i][0] = DD.x[i];
			DDCopy[i][1] = DD.y[i];
			DDCopy[i][2] = DD.z[i];
			DDCopy[i][3] = ((Rb + 0.0)*(Rb + 0.0))*DD.q[i] * Sigma_rr[i] / (4 * pi);
	                  DDDCopy[i][0] = DDD.x[i];
			DDDCopy[i][1] = DDD.y[i];
			DDDCopy[i][2] = DDD.z[i];
			DDDCopy[i][3] = ((Rb + 0.0)*(Rb + 0.0))*DDD.q[i] * Sigma_rr[i] / (4 * pi);
		} 
		//GPUset
		double4 *SourceD,*SourceDD,*SourceDDD;
		cudaMalloc((void**)&SourceD,sizeof(double4)*Np);
		cudaMemcpy(SourceD,DCopy,sizeof(double4)*Np,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&SourceDD,sizeof(double4)*Np);
		cudaMemcpy(SourceDD,DDCopy,sizeof(double4)*Np,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&SourceDDD,sizeof(double4)*Np);
		cudaMemcpy(SourceDDD,DDDCopy,sizeof(double4)*Np,cudaMemcpyHostToDevice);
		
		number[0] = N;
		number[1] = Np;
		cudaMemcpy(Number1, number, sizeof(int) * 2, cudaMemcpyHostToDevice);

		BlockNum = N / 128 + 1;
		ThreadNum = 128;
 
    double *GPUFieldPotentialD,*GPUFieldPotentialDD,*GPUFieldPotentialDDD;
		cudaMalloc((void**)&GPUFieldPotentialD, sizeof(double)*N);
    cudaMalloc((void**)&GPUFieldPotentialDD, sizeof(double)*N);
    cudaMalloc((void**)&GPUFieldPotentialDDD, sizeof(double)*N);

		//Calculate
		calculate_forces1 << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (PositionField1, SourceD, GPUFieldPotentialD, Number1);
		calculate_forces1 << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (PositionField1, SourceDD, GPUFieldPotentialDD, Number1);
    calculate_forces1 << <BlockNum, ThreadNum, sizeof(double4)*ThreadNum >> > (PositionField1, SourceDDD, GPUFieldPotentialDDD, Number1);	
			                 
		//Copy from GPU to CPU
    //near field potential
    cudaStatus = cudaSetDevice(0); 

    cudaMemcpy(&GPUFIELDPOTENTIALSECOND, GPUFieldPotentialSecond, sizeof(double)*qq, cudaMemcpyDeviceToHost);
		cudaFree(SourceSecond);
		cudaFree(PositionFieldSecond);
		cudaFree(NumberSecond);
		cudaFree(GPUFieldPotentialSecond);
    cudaStatus = cudaSetDevice(1); 
		cudaMemcpy(&GPUFIELDPOTENTIAL, GPUFieldPotential, sizeof(double)*qq, cudaMemcpyDeviceToHost);
		cudaFree(Source);
		cudaFree(PositionField);
		cudaFree(Number);
    cudaFree(GPUFieldPotential);  
        
		finish1 = clock();
    double duration1;
		duration1 = (double)(finish1 - start1+0.00) / (CLOCKS_PER_SEC+0.00);
		printf("Time for GPU to calculate the near field potential is:");
		printf("%f seconds\n", duration1);
	                 
		double GPUFIELDPOTENTIALFARD[N],GPUFIELDPOTENTIALFARDD[N],GPUFIELDPOTENTIALFARDDD[N];
		cudaMemcpy(&GPUFIELDPOTENTIALFARD, GPUFieldPotentialD, sizeof(double)*N, cudaMemcpyDeviceToHost);
    cudaMemcpy(&GPUFIELDPOTENTIALFARDD, GPUFieldPotentialDD, sizeof(double)*N, cudaMemcpyDeviceToHost); 
    cudaMemcpy(&GPUFIELDPOTENTIALFARDDD, GPUFieldPotentialDDD, sizeof(double)*N, cudaMemcpyDeviceToHost); 
                   
		cudaFree(SourceD);
    cudaFree(SourceDD);
    cudaFree(SourceDDD);

		cudaFree(Number1);
		cudaFree(PositionField1);

		finish1 = clock();
		printf("Time for GPU to calculate the far field potential is:");
		duration = ((double)(finish1 - start1+0.00)) / (CLOCKS_PER_SEC+0.00)-duration;
		printf("%f seconds\n", duration);  
           
		/*       end   CUDA            */
  
		start1=clock();

		double tt = 0.00;             
		for (int i = 0; i < N; i++)
		{
			
			tt = tt +(GPUFIELDPOTENTIALFARD[i]- GPUFIELDPOTENTIALFARDD[i]+GPUFIELDPOTENTIALFARDDD[i]);
		}
		for(int i=0;i<qq;i++)
		{
		  tt=tt+GPUFIELDPOTENTIAL[i]+GPUFIELDPOTENTIALSECOND[i]; 
		}
                   
                  
		cout << "Calculated by simulate is " << setprecision(16) << tt / (2.0) << endl;    //This is the result for Madelung
		cout << "The Absolute Error is" << abs(-470125.129478090000 - tt / 2.0) << endl;
		cout << "The Relative Error is" << abs(1 - (-470125.129478090000 / (tt / 2.0))) << endl;

    finish1 = clock();
		printf("Time for GPU to calculate the total sum is:");
		duration = (double)(finish1 - start1) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration);                

		finishFULL = clock();
		printf("Time for all of this prodecure is:");
		duration = (double)(finishFULL - startFULL) / CLOCKS_PER_SEC;
		printf("%f seconds\n", duration-cut1);
                                       
  	return 1;
}

