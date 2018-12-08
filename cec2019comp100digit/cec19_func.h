#ifndef _CEC19_FUNC_H_
#define _CEC19_FUNC_H_

void cec19_test_func(double *, double *,int,int,int);

#ifdef _MAIN
#define EXTERN extern
#else
#define EXTERN 
#endif

EXTERN double *OShift,*M,*y,*z,*x_bound;
EXTERN int ini_flag,n_flag,func_flag,*SS;


#endif
