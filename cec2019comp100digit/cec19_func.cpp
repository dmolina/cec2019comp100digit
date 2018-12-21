/*
  CEC19 Test Function Suite 
  Noor Awad (email: noorawad1989@gmail.com) 
  Dec. 13th 2018
  1. Run the following command in Matlab window:
  mex cec19_func.cpp -DWINDOWS
  2. Then you can use the test functions as the following example:
  f = cec19_func(x,func_num); 
  Here x is a D*pop_size matrix.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "cec19_func.h"

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void Lennard_Jones(double *, int, double *); /* Lennard Jones */
void Hilbert(double *, int, double *); /* Hilbert */
void Chebyshev(double *, int, double *); /* Chebyshev */
void schaffer_F7_func (double *, double *, int , double *,double *, int, int); /* Schwefel's F7 */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void weierstrass_func (double *, double *, int , double *,double *, int, int); /* Weierstrass's  */
void schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer¡¯s F6  */
void happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */

void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void asyfunc (double *, double *x, int, double);
void oszfunc (double *, double *, int);

void cec19_test_func(double *, double *,int,int,int);

static const char *dir_path;

void set_dir_path(const char *path) {
  dir_path = strdup(path);
}

void cec19_test_func(double *x, double *f, int nx, int mx,int func_num)
{
	int i;
	if (ini_flag==1)
	{
    assert(dir_path != nullptr);

		if ((n_flag!=nx)||(func_flag!=func_num))
		{
			ini_flag=0;
		}
	}

	if (ini_flag==0)
	{
		FILE *fpt;
		char FileName[256];
		free(M);
		free(OShift);
		free(y);
		free(z);
		free(x_bound);
		y=(double *)malloc(sizeof(double)  *  nx);
		z=(double *)malloc(sizeof(double)  *  nx);
		x_bound=(double *)malloc(sizeof(double)  *  nx);
		for (i=0; i<nx; i++)
			x_bound[i]=100.0;

		if (!(nx==2||nx==10||nx==9||nx==16||nx==18))
		{
			printf("Error: Test functions are only defined for D=10, 9, 16, 18 \n F1 is defined on D=9 \n F2 is defined on D=16 \n F3 is defined on D=18 \n F4-F10 are defined on D=10.");
		}
		

		/* Load Matrix M*/
		if (func_num > 3)
		{
      sprintf(FileName, "%s/M_%d_D%d.txt", dir_path, func_num,nx);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
      printf("\n Error: Cannot open input '%s' file for reading \n", FileName);
		}
		M=(double*)malloc(nx*nx*sizeof(double));
		if (M==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for (i=0; i<nx*nx; i++)
		{
			fscanf(fpt,"%lf",&M[i]);
		}
		
		fclose(fpt);
		}
		
		/* Load shift_data */
		if (func_num > 3)
		{
      sprintf(FileName, "%s/shift_data_%d.txt", dir_path, func_num);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input '%s' file for reading \n", FileName);
		}

		OShift=(double *)malloc(nx*sizeof(double));
		if (OShift==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for(i=0;i<nx;i++)
		{
			fscanf(fpt,"%lf",&OShift[i]);
		}
		
		fclose(fpt);
		}
		
		n_flag=nx;
		func_flag=func_num;
		ini_flag=1;
		//printf("Function has been initialized!\n");
	}


	for (i = 0; i < mx; i++)
	{
		switch(func_num)
		{
		case 1:
			Chebyshev(&x[i*nx], nx, &f[i]);	
			f[i]+=1.0;
			break;
			
		case 2:	
			Hilbert(&x[i*nx], nx, &f[i]);	
			f[i]+=1.0;
			break;
			
		case 3:	
			Lennard_Jones(&x[i*nx], nx, &f[i]);	
			f[i]+=1.0;
			break;
			
		case 4:	
			rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 5:	
			griewank_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 6:	
			weierstrass_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 7:	
			schwefel_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 8:
			escaffer6_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 9:
			happycat_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
			
		case 10:	
			ackley_func(&x[i*nx],&f[i],nx,OShift,M,1,1);
			f[i]+=1.0;
			break;
		default:
			printf("\nError: There are only 30 test functions in this test suite!\n");
			f[i] = 0.0;
			break;
		}
		
	}

}


void schaffer_F7_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's 1.2  */
{
    int i;
	double tmp;
    f[0] = 0.0;
	sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
	for (i=0; i<nx-1; i++)	
	{
		z[i]=pow(y[i]*y[i]+y[i+1]*y[i+1],0.5);
		tmp=sin(50.0*pow(z[i],0.2));
		f[0] += pow(z[i],0.5)+pow(z[i],0.5)*tmp*tmp ;
	}
	f[0] = f[0]*f[0]/(nx-1)/(nx-1);
}

void griewank_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;

	sr_func (x, z, nx, Os, Mr, 600.0/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		s += z[i]*z[i];
		p *= cos(z[i]/sqrt(1.0+i));
	}
	f[0] = 1.0 + s/4000.0 - p;
}

void ackley_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;

	sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		sum1 += z[i]*z[i];
		sum2 += cos(2.0*PI*z[i]);
	}
	sum1 = -0.2*sqrt(sum1/nx);
	sum2 /= nx;
		f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}

void weierstrass_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    f[0] = 0.0;

	sr_func (x, z, nx, Os, Mr, 0.5/100.0, s_flag, r_flag); /* shift and rotate */

  assert(nx > 0);

	for (i=0; i<nx; i++)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (j=0; j<=k_max; j++)
		{
			sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(z[i]+0.5));
			sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
		}
		f[0] += sum;
	}
	f[0] -= nx*sum2;
}


void rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rastrigin's  */
{
    int i;
	f[0] = 0.0;

	sr_func (x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		f[0] += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
	}
}

void step_rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
{
    int i;
	f[0]=0.0;
	for (i=0; i<nx; i++)
	{
		if (fabs(y[i]-Os[i])>0.5)
		y[i]=Os[i]+floor(2*(y[i]-Os[i])+0.5)/2;
	}

	sr_func (x, z, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		f[0] += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
	}
}

void schwefel_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's  */
{
    int i;
	double tmp;
	f[0]=0.0;

	sr_func (x, z, nx, Os, Mr, 1000.0/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		z[i] += 4.209687462275036e+002;
		if (z[i]>500)
		{
			f[0]-=(500.0-fmod(z[i],500))*sin(pow(500.0-fmod(z[i],500),0.5));
			tmp=(z[i]-500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else if (z[i]<-500)
		{
			f[0]-=(-500.0+fmod(fabs(z[i]),500))*sin(pow(500.0-fmod(fabs(z[i]),500),0.5));
			tmp=(z[i]+500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else
			f[0]-=z[i]*sin(pow(fabs(z[i]),0.5));
		}
		f[0] +=4.189828872724338e+002*nx;

}



void escaffer6_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Expanded Scaffer¡¯s F6  */
{
    int i;
    double temp1, temp2;

	sr_func (x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

    f[0] = 0.0;
    for (i=0; i<nx-1; i++)
    {
        temp1 = sin(sqrt(z[i]*z[i]+z[i+1]*z[i+1]));
		temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z[i]*z[i]+z[i+1]*z[i+1]);
        f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z[nx-1]*z[nx-1]+z[0]*z[0]));
	temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z[nx-1]*z[nx-1]+z[0]*z[0]);
    f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
}

void happycat_func (double *x, double *f, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
{
	int i;
	double alpha,r2,sum_z;
	alpha=1.0/8.0;
	
	sr_func (x, z, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */

	r2 = 0.0;
	sum_z=0.0;
    for (i=0; i<nx; i++)
    {
		z[i]=z[i]-1.0;//shift to orgin
        r2 += z[i]*z[i];
		sum_z += z[i];
    }
    f[0]=pow(fabs(r2-nx),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
}


void shiftfunc (double *x, double *xshift, int nx,double *Os)
{
	int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc (double *x, double *xrot, int nx,double *Mr)
{
	int i,j;
    for (i=0; i<nx; i++)
    {
        xrot[i]=0;
			for (j=0; j<nx; j++)
			{
				xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
			}
    }
}

void sr_func (double *x, double *sr_x, int nx, double *Os,double *Mr, double sh_rate, int s_flag,int r_flag) /* shift and rotate */
{
	int i;
	if (s_flag==1)
	{
		if (r_flag==1)
		{	
			shiftfunc(x, y, nx, Os);
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				y[i]=y[i]*sh_rate;
			}
			rotatefunc(y, sr_x, nx, Mr);
		}
		else
		{
			shiftfunc(x, sr_x, nx, Os);
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				sr_x[i]=sr_x[i]*sh_rate;
			}
		}
	}
	else
	{	

		if (r_flag==1)
		{	
			for (i=0; i<nx; i++)//shrink to the orginal search range
			{
				y[i]=x[i]*sh_rate;
			}
			rotatefunc(y, sr_x, nx, Mr);
		}
		else
		for (i=0; i<nx; i++)//shrink to the orginal search range
		{
			sr_x[i]=x[i]*sh_rate;
		}
	}
}

void asyfunc (double *x, double *xasy, int nx, double beta)
{
	int i;
    for (i=0; i<nx; i++)
    {
		if (x[i]>0)
        xasy[i]=pow(x[i],1.0+beta*i/(nx-1)*pow(x[i],0.5));
    }
}

void oszfunc (double *x, double *xosz, int nx)
{
  const double INIT_XX = -1000;
	int i,sx;
	double c1,c2,xx=INIT_XX;

    for (i=0; i<nx; i++)
    {
		if (i==0||i==nx-1)
        {
			if (x[i]!=0)
				xx=log(fabs(x[i]));
			if (x[i]>0)
			{	
				c1=10;
				c2=7.9;
			}
			else
			{
				c1=5.5;
				c2=3.1;
			}	
			if (x[i]>0)
				sx=1;
			else if (x[i]==0)
				sx=0;
			else
				sx=-1;
      assert(xx != INIT_XX);
			xosz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
		}
		else
			xosz[i]=x[i];
    }
}



void Lennard_Jones(double *x,int D, double *f)  // find the atomic configuration with minimum energy
{
	/* valid for any dimension, D=3*k, k=2,3,4,...,25.   k is the number of atoms in 3-D space
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 4, initial lower bound = -4
	value-to-reach = minima[k-2]+.0001
	f(x*) = minima[k-2]; see array of minima below; additional minima available at the
	Cambridge cluster database: http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html
	*/

	f[0] = 0;
	int i, j, k, a, b;
	long double xd, yd, zd, ed, ud, sum = 0;


	k = D / 3;
	if (k < 2)  // default if k<2
	{
		k = 2;
		D = 6;
	}

	for (i = 0; i < k - 1; i++)
	{
		for (j = i + 1; j < k; j++)
		{
			a = 3 * i;
			b = 3 * j;
			xd = x[a] - x[b];
			yd = x[a + 1] - x[b + 1];
			zd = x[a + 2] - x[b + 2];
			ed = xd*xd + yd*yd + zd*zd;
			ud = ed*ed*ed;
			if (ud > 1.0e-10) sum += (1.0 / ud - 2.0) / ud;
			else sum += 1.0e20; 
		}
	}
	
	f[0] += sum;
	f[0] += 12.7120622568;
}

void Hilbert(double *x, int D, double *f)  // find the inverse of the (ill-conditioned) Hilbert matrix
{
	/* valid for any dimension, n=k*k, k=2,3,4,...
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^n, initial lower bound = -2^n
	value-to-reach = f(x*)+1.0e-8
	f(x*) = 0.0; x*={{9,-36,30},{-36,192,-180},{30,-180,180}} (n=9)
	x*={{16,-120,240,-140},{-120,1200,-2700,1680},{240,-2700,6480,4200},{-140,1680,-4200,2800}} (n=16)
	*/

	f[0] = 0;
	int i, j, k, b;

	long double sum = 0;

	static long double hilbert[10][10], y[10][10];			// Increase matrix size if D > 100

	b = (int)sqrt((double)D);

	
	for (i = 0; i < b; i++)
		{
			for (j = 0; j < b; j++)
			{
				hilbert[i][j] = 1. / (double)(i + j + 1);		// Create a static Hilbert matrix
			}
		}

	for (j = 0; j < b; j++)
	{
		for (k = 0; k < b; k++)
		{
			y[j][k] = 0;
			for (i = 0; i < b; i++)
			{
				y[j][k] += hilbert[j][i] * x[k + b * i];		// Compute matrix product H*x
			}
		}
	}


	for (i = 0; i < b; i++)
	{
		for (j = 0; j < b; j++)
		{
			if (i == j) sum += fabs(y[i][j] - 1);				// Sum absolute value of deviations
			else sum += fabs(y[i][j]);
		}
	}
	
	f[0] += sum;
}

void Chebyshev( double *x, int D, double *f)  // Storn's Tchebychev - a 2nd ICEO function - generalized version
{
	/* Valid for any D>2
	constraints: unconstrained
	type: multi-modal with one global minimum; non-separable
	initial upper bound = 2^D, initial lower bound = -D^n
	value-to-reach = f(x*)+1.0e-8
	f(x*)=0.0; x*=(128,0,-256,0,160,0,-32,0,1) (n=9)
	x*=(32768,0,-131072,0,212992,0,-180224,0,84480,0,-21504,0,2688,0,-128,0,1) (n=17)
	*/

	f[0] = 0.0;
	int i, j;
	static int sample;
	long double a = 1., b = 1.2, px, y = -1, sum = 0;
	static long double dx, dy;
	
	for (j = 0; j < D - 2; j++)
		{
			dx = 2.4 * b - a;
			a = b;
			b = dx;
		}

		sample = 32 * D;
		dy = 2. / (long double)sample;

	for (i = 0; i <= sample; i++)
	{
		px = x[0];
		for (j = 1; j < D; j++)
		{
			px = y*px + x[j];
		}
		if (px < -1 || px > 1) sum += (1. - fabs(px))*(1. - fabs(px));
		y += dy;
	}

	for (i = -1; i <= 1; i += 2)
	{
		px = x[0];
		for (j = 1; j < D; j++)
		{
			px = 1.2*px + x[j];
		}

		if (px < dx) sum += px * px;
	}
	
	f[0] += sum;
	
}



