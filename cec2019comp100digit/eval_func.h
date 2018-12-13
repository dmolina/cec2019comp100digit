#ifndef _CEC19_EVAL_FUN_H

#define _CEC19_EVAL_FUN_H

void set_func(int funcID, int dim);
void set_dir_path(const char *new_data_dir);
double eval_sol(double *x);
void free_func(void);

#endif
