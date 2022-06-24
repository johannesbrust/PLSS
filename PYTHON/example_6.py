"""
Test PLSS (Underdetermined)

Projected Linear Systems Solver
"""
#--------------------------- example_6.py -------------------------------------
#
# Test of the Python 3 versions of PLSS on an underdetermined system
#
#------------------------------------------------------------------------------
# J.B., 05/10/22, Initial version
# J.B., 05/22/22, Test of Kaczmarz like iteration
# J.B., 06/24/22, Underdetermined test

# Import
import numpy as np
import plss_KZ as plkz
import plss_R as plr
import plss_RW2 as plrw

# Test uses a matrix from SOL (Systems Optimization Laboratory)
# in testing CRAIG's method

n = 6;
m = 3;

# Form the "kahan" matrix (upper trapezoidal)
A       = np.array(
        [[1.0000,   -0.3624,   -0.3624,   -0.3624,   -0.3624,   -0.3624],
         [0,    0.9320,   -0.3377,   -0.3377,   -0.3377,   -0.3377],
         [0,         0,    0.8687,   -0.3148,   -0.3148,   -0.3148]]);
xsol    = np.ones(n);

b = A @ xsol

# Call PLSS KZ
# All default values are used, except results are printed
x0 = np.zeros(np.size(xsol))
(xk,rk,outs) = plkz.plss_KZ(x0,A,b,maxiter=n+1,printF=1)

# Call PLSS R
# All default values are used, except results are printed
(xkr,rkr,outsr) = plr.plss_R(x0,A,b,maxiter=n+1,printF=1)

# Call PLSS RW
# All default values are used, except results are printed
(xkrw,rkrw,outsrw) = plrw.plss_RW2(x0,A,b,maxiter=n+1,printF=1)