#!python
#cython: language_level=2, boundscheck=False
from os import path
from collections import namedtuple
from pkg_resources import resource_filename
import cython

cdef extern from "eval_func.h":
    void set_func(int funid, int dim)
    double eval_sol(double*)
    void set_dir_path(const char * new_data_dir)
    void free_func()

def init(int fun, int dim):
    """
    Evaluate the solution
    """
    free_func()
    set_func(fun, dim)
    cdef bytes dir_name = resource_filename(__name__, "input_data").encode()
    set_dir_path(dir_name)

def eval(double[::1] x):
    fitness = eval_sol(&x[0])
    return fitness

def end():
    free_func()
