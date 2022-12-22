import numpy as np
import time as time

def plss_RSP(x, A, b,
 tol=1e-6,
 maxiter=500,
 hasMeasure=0,
 measure=0,
 printF=0,
 store=0):
    """Projected linear systems solver"""
    #PLSS_RSP: Projected linear systems solver Algorithm 1 (Recursive)
    # Sparse computations
    # PYTHON 3.8 Implementation
    # Method for iteratively solving linear systems
    # 
    # Ax = b,
    #
    # where A in (mxn), x in (nx1), b in (mx1).
    #
    # This algorithm computes search directions pk defined by the projections
    # defined by sketches that store previous residuals in -ck = rk = b - A*xk.
    # 
    # Let Sk = [c0,...,ck], then the step is
    # 
    # pk = -A Sk (Sk' AA Sk)^{-1} (Sk' Sk) e_k,                             (1)
    # 
    # where e_k denotes the kth column of the identity matrix. 
    # When (1) is simplified, it can be expressed by the recursive formula
    #
    # pk = beta1*p_{k-1} + beta_2*yk,                                       (2)
    #
    # where yk = A'ck.
    #
    # INPUTS:
    # x := Initial guess (nx1)
    # A := Linear system (mxn)
    # b := Right hand side (mx1)
    # tol := Residual tolerance 
    # maxiter  := Maximum iterations
    # hasMeasure:= Flag for error measure (for comparisons). 
    # measure:= Different error measure (for comparisons). 
    # printF  := Flag for printing
    # store  := Flag for storing information (for comparisons)
    #
    # OUTPUTS:
    # xk := Iterate at termination
    # nck:= norm(ck=(A*xk-b)), 
    # outs:= Dictionary with solver outpus
    #   outs['ex'] := Flag is converged to tolerance (1 conv, 0 not)
    #   outs['ctime']:= Compute time used by solver
    #   outs['niter']:= Number of iterations
    #   outs['numA']:= Number of matrix multiplies A*(VEC)
    #   [outs['errs']]:= [Optional] If opts.store==true, then errors stored
    #   [outs['times']]:= [Optional] If opts.store==true, then times stored
    #-------------------------------------------------------------------------#
    # 04/02/20, J.B., Initial implementation
    # 12/07/20, J.B., Recursive formula
    # 05/08/22, J.B., Preparation for release
    # 05/09/22, J.B., Python implementation
    # 12/22/22, J.B., Sparse computation
    
    # Initializing storage    
    tstart  = time.time();
    
    #A = aslinearoperator(A);
    
    numA    = 0;
    ex      = 0;
    
    if store == 1:      
        
        errs  = [];
        times = [];
    
    xk      = x;
    ck      = A.dot(xk)-b;
    numA    = numA + 1;
    
    nck     = np.linalg.norm(ck);
    yk      = A.transpose().dot(ck/nck);
    
    numA    = numA + 1;
    k       = 0;
    
    # Store minimum solution estimate
    xkmin = xk;
    nckmin = nck;
    
    if hasMeasure == 1:
        ckm = np.sum(measure*ck,1);
        nck = np.sqrt(np.sum((ckm*ckm)));
    else:
        nck = np.linalg.norm(ck);
    
    if store == 1:
        errs.append(nck);
        times.append(time.time()-tstart);
    
    
    if nck <= tol:                
        ex      = 1;
        tend    = time.time()-tstart;
        
    
    rhok        = nck;
    deltaki     = 1/np.sum((yk*yk));
    pk          = -(deltaki*rhok)*yk;
    
    if printF == 1:
    
        (m,n) = A.shape;
        
        print('');
        print('********* Linear Systems Solver ********************  ');
        print('*         Alg: PLSS (Python 3.9)                      ');
        print('*         Size: m = %i, n = %i                        '% (m,n));
        print('*         tol = %1.8f                                 '% tol);
        print('*         Maxit = %i                                  '% maxiter);        
        print('****************************************************  ');        
        print('');
        
        print('Iter \t norm(Res) \t norm(pk)    ');
        print('%i \t %.4e \t %.4e   ' % (k,nck,np.linalg.norm(pk)));
        
   
    
    k           = k + 1;
    xk          = xk + pk;
    
    while np.logical_and(np.logical_and(tol < nck,k < maxiter),ex == 0):
            
        ck      = A.dot(xk) - b;     numA = numA + 1;    
        
        nck     = np.linalg.norm(ck);
        
        yk      = A.transpose().dot(ck/nck);  numA = numA + 1;
            
        rhok    = nck;    
        
        p2      = np.sum((pk*pk));
        nrp      = np.sqrt(p2);    
        
        py      = np.sum((pk*yk));
        yy      = np.sum((yk*yk));
        ny      = np.sqrt(yy);
        
        denom   = (nrp*ny-py)*(nrp*ny+py);
        
        beta1   = (rhok * py)/denom;
        beta2   = -(rhok * p2)/denom;
        
        # Step computation
        pk      = beta1*pk + beta2*yk;
            
        # Iteration information
        if hasMeasure == 1:
            ckm = np.sum(measure*ck,1);
            nck = np.sqrt(np.sum((ckm*ckm)));
        else:
            nck = np.linalg.norm(ck); 
            
        # Safeguard (to set pk=0 when a solution is found)
        p2      = np.sum((pk*pk));
        nrp     = np.sqrt(p2);
        if np.logical_and(nck <= tol,np.isnan(nrp)):
            pk[:] = 0;
            nrp = 0;
                    
        if np.logical_and(printF == 1,np.logical_or(k < 11,np.mod(k,20)==0)):
            print('%i \t %.4e \t %.4e  ' % (k,nck,nrp));         
        

        if store == 1:
            errs.append(nck);
            times.append(time.time()-tstart);
        
        
        if nck <= tol:        
            ex          = 1; 
            tend        = time.time()-tstart;            
            if np.logical_and(np.logical_and(printF == 1,k > 11),np.mod(k,20)!=0):
                print('%i \t %.4e \t %.4e  ' % (k,nck,np.sqrt(np.sum((pk*pk)))));
                               
        # Store minimum iterate
        if nckmin >= nck:            
            xkmin = xk;
            nckmin = nck;
           
        # Prepare for next iteration
        xk  = xk + pk;
        k   = k + 1;
    
    
    if ex == 0:
        tend = time.time()-tstart;    
    
    rk = A.dot(xk)-b;
    nck = np.sqrt(np.sum((rk*rk)));
    
    numA = numA+1;
    
    if nck < tol:
        ex = 1;

    
    if nckmin < nck:
       
        nck = nckmin;
        xk = xkmin;
    
    outs = {'ex': ex,'ctime': tend, 'niter': k, 'numA': numA};
        
    if store == 1:        
        outs['errs']   = errs;
        outs['times']  = times;
    
    if printF == 1:
        
        print('');
        print('********* Summary **********************************  ');
        print('*         Conv: %i                                    ' % ex);
        print('*         Time (s): %1.2f                             ' % tend);    
        print('*         norm(Res) = %1.4f                           ' % nck);
        print('*         Iter = %i                                   ' % k);
        print('*         Num. mult. = %i                             ' % numA);
        print('****************************************************  ');
        print('');
                
    return  xk,nck,outs
