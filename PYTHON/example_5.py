"""
Test PLSS (Kaczmarz Method)

Projected Linear Systems Solver
"""
#--------------------------- example_5.py -------------------------------------
#
# Test of the Python 3 version of PLSS KZ
#
# This problem is based on transposing the bidiagonal matrix from 
# experiments using CRAIG's method from Stanford SOL
# 
# The solution with n=10 is [n,n-1,...,1]
#------------------------------------------------------------------------------
# J.B., 05/10/22, Initial version
# J.B., 05/22/22, Test of Kaczmarz like iteration

# Import
import numpy as np
import plss_KZ as plkz

# Test uses a matrix from SOL (Systems Optimization Laboratory)
# in testing CRAIG's method

n = 10;
m = n;

# Form bidiagonal matrix and solution
A       = np.diag(np.arange(1., 1.+n)) + np.diag(np.arange(1., n),-1);
xsol    = np.array(np.arange(n,0,-1.));

b = A.T @ xsol

# Call PLSS KZ
# All default values are used, except results are printed
x0 = np.zeros(np.size(xsol))

(xk,rk,outs) = plkz.plss_KZ(x0,A.T,b,maxiter=n+1,printF=1)


