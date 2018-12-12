#define _MAIN
#include "cec19_func.h"
#include <stdlib.h>

static int _funcID;
static int _dim;

void set_func(int funcID, int dim) {
  _funcID = funcID;
  _dim = dim;
  ini_flag = 0;
}

double eval_sol(double *x) {
  double fit;
  cec19_test_func(x, &fit, _dim, 1, _funcID);
  return fit;
}

void free_func(void) {
	free(z);
	free(M);
	free(OShift);
	free(x_bound);
  // Init to allow set_func again
  z = nullptr;
  M = nullptr;
  OShift = nullptr;
  x_bound = nullptr;
  ini_flag = 0;
}
