"""
Test PLSS

Projected Linear Systems Solver
"""
#--------------------------- example_7.py ------------------------------------
#
# Test of the Python 3 version of PLSS
#
# This problem is a bidiagonal matrix from experiments using 
# CRAIG's method from Stanford SOL
# 
# The solution with n=10 is [n,n-1,...,1]
#
# This example is the result of a question, based an a performance
# comparison raised on  Dec. 22, 2022 on github (asked by Vandermode)
#
# Based on the question sparse versions of plss_R and plss_RW2 were added
# These two versions can be tried in this experiment
#
#------------------------------------------------------------------------------
# J.B., 12/22/22, Initial version

import numpy as np
import plss_RSP as pl
#import plss_RW2SP as pl

from scipy.sparse import csr_matrix

# Test uses a matrix from SOL (Systems Optimization Laboratory)
# in testing CRAIG's method

n = 1000
m = n

# Form bidiagonal matrix and solution
A       = np.diag(np.arange(1., 1.+n)) + np.diag(np.arange(1., n),-1);
xsol    = np.array(np.arange(n,0,-1.));

ASP     = csr_matrix(A);

b = A.dot(xsol);

# Call PLSS
# All default values are used, except results are printed
x0 = np.zeros(np.size(xsol))

xk_np = np.linalg.solve(A, b)
print(xk_np)
(xk,rk,outs) = pl.plss_RSP(x0,ASP,b,printF=0, maxiter=100000)
#(xk,rk,outs) = pl.plss_RW2SP(x0,ASP,b,printF=0, maxiter=100000)
print(xk)
print(rk)
print(outs)