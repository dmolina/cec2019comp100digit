/*
  CEC19 Test Function Suite for Single Objective Optimization
  Noor Awad (email: noorawad1989@gmail.com) 
  Sep. 21th 2018
*/

#include <stdio.h>
#include <math.h>
#include <malloc.h>


void cec17_test_func(double *, double *,int,int,int);

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;


int main(void)
{
	int i,j,k,n,m,func_num;
	double *f,*x;
	FILE *fpt;
	char FileName[10];
	m=2;
	n=10;
	x=(double *)malloc(m*n*sizeof(double));
	f=(double *)malloc(sizeof(double)  *  m);
	for (i = 0; i < 10; i++)
	{
		func_num=i+1;
		sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		
		if (x==NULL)
			printf("\nError: there is insufficient memory available!\n");

		for(k=0;k<n;k++)
		{
				fscanf(fpt,"%lf",&x[k]);
				/*printf("%Lf\n",x[k]);*/
		}

		fclose(fpt);

			for (j = 0; j < n; j++)
			{
				x[1*n+j]=0.0;
				/*printf("%Lf\n",x[1*n+j]);*/
			}
		
		
		for (k = 0; k < 1; k++)
		{
			cec19_test_func(x, f, n,m,func_num);
			for (j = 0; j < 2; j++)
			{
				printf(" f%d(x[%d]) = %lf,",func_num,j+1,f[j]);
			}
			printf("\n");
		}
	
	}
	free(x);
	free(f);
	free(y);
	free(z);
	free(M);
	free(OShift);
	free(x_bound);
  return 0;
}


