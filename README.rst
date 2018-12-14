Introduction
============

This is a Python wrapping using the C++ Implementation of the test suite for the
Special Session on Large Scale Global Optimization at 2019 IEEE Congress on
Evolutionary Computation http://cec2019.org/programs/competitions.html#cec-06.

http://www.ntu.edu.sg/home/epnsugan/index_files/CEC2019/CEC2019.htm


Note
----
If you are to use any part of this code, please cite the following publications:

   K. V. Price, N. H. Awad, M. Z. Ali, P. N. Suganthan, "Problem Definitions and
   Evaluation Criteria for the 100-Digit Challenge Special Session and
   Competition on Single Objective Numerical Optimization,"  Technical Report,
   Nanyang Technological University, Singapore, November 2018.

http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%2d2019&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

Requirements
------------

- GNU Make
- GNU G++
- Python
- Cython

Testing Environment
-------------------

- Debian GNU/Linux jessie/sid
- GNU Make 3.81
- g++ (Debian 4.7.3-4) 4.7.3
- Python 2.7 and Python 3.2
- numpy 1.8.1
- cython 0.20.1

Instalation
-----------

It is pending to submit to pip, when it is ready.

Very easy, *pip install cec2019comp100digit* ;-). 

You can also download from https://github.com/dmolina/cec2019comp100digit, and do *python setup.py install [--user]*.
(the option *--user* is for installing the package locally, as a normal user (interesting when you want to 
run the experiments in a cluster/server without administration permissions).

To compile the source code in C++
----------------------------------

The source code in C++ is also available. If you want to compile only the C++
version type in 'make' in the root directory of source code. 

There are two equivalents demo executables: demo and demo2. 

**REMEMBER: To run the C++ version the directory input_data must be available in the working directory**. 
In the python version, these files are included in the packages, so it is not
needed.

Quickstart
----------

The package is very simple to use. There is a package cec2019comp100digit with
three functions:

- **init(fun_id, Dim)**
  Init the function for the dimension selected.

- **eval(sol)**
  Eval the solution, when sol is a numpy (or array) of dimension *Dim*.

- **end()**
  Free resources.

Init function
-------------
>>> from cec2019comp100digit import cec2019comp100digit
>>> bench = cec2019comp100digit
>>> bench.init(3, 10) # Init function 3

Create a random solution
~~~~~~~~~~~~~~~~~~~~~~~~
>>> import numpy as np
>>> sol = np.random.rand(10)

Evaluate a solution
~~~~~~~~~~~~~~~~~~~
>>> bench.eval(sol)
18010038.104525752

Freeing resources
~~~~~~~~~~~~~~~~~
>>> bench.end()

Contact
-------

Python package 
  Daniel Molina @ Computer Science Deparment, University of Granada
  Please feel free to contact me at <dmolina@decsai.ugr.es> for any enquiries or
  suggestions.


Last Updated
~~~~~~~~~~~~
- C++ version
  <2018-12-08>

- Python wrapping
  <2018-12-08>
