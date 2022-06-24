"""
Test PLSS

Projected Linear Systems Solver
"""
#--------------------------- example_3.py -------------------------------------
#
# Test of the Python 3 version of PLSS
#
# This problem is a bidiagonal matrix from experiments using 
# CRAIG's method from Stanford SOL
# 
# The solution with n=10 is [n,n-1,...,1]
#------------------------------------------------------------------------------
# J.B., 05/10/22, Initial version

# Import
import numpy as np
import plss_R as pl

# Test uses a matrix from SOL (Systems Optimization Laboratory)
# in testing CRAIG's method

n = 10;
m = n;

# Form bidiagonal matrix and solution
A       = np.diag(np.arange(1., 1.+n)) + np.diag(np.arange(1., n),-1);
xsol    = np.array(np.arange(n,0,-1.));

b = A @ xsol

# Call PLSS
# All default values are used, except results are printed
x0 = np.zeros(np.size(xsol))

(xk,rk,outs) = pl.plss_R(x0,A,b,printF=1)