import numpy as np
import time as time

def plss_KZ(x, A, b,
 tol=1e-6,
 maxiter=500,
 hasMeasure=0,
 measure=0,
 printF=0,
 store=0):
    """Projected linear systems solver Kaczmarz"""
    #PLSS_KZ_EX: Projected linear systems solver using a Kaczmark like process. 
    #
    # Method for iteratively solving linear systems
    # 
    # Ax = b,
    #
    # where A in (mxn), x in (nx1), b in (mx1).
    #
    # This algorithm computes search directions pk defined by the projections
    # defined by sketches that store identity columns.
    # 
    # Let Sk = [e1,...,ek1], then the step is
    # 
    # pk = -A Sk (Sk' AA Sk)^{-1} (Sk' ck),
    # 
    # where e_k denotes the kth column of the identity matrix. This method
    # stores the updates Pk=[p0,...,pkm1] and the products PkA =
    # [A*p0,...,A*pkm1]. 
    # Options for permutation of rows is possible:
    #       perm =  0 : No reordering
    #               1 : Increasing row norms
    #               2 : Decreasing residuals
    #               3 : dmperm "Dulmage-Mendelsohn permutation"
    #               4 : randomized (default)
    # Note, on restart the iteration starts with decreasing residuals
    #
    # INPUTS:
    # x := Initial guess (nx1)
    # A := Linear system (mxn)
    # b := Right hand side (mx1)
    # maxiter  := Maximum iterations
    # hasMeasure:= Flag for error measure (for comparisons). 
    # measure:= Different error measure (for comparisons). 
    # printF  := Flag for printing
    # store  := Flag for storing information (for comparisons)
    #
    # OUTPUTS:
    # xk := Iterate at termination
    # nck:= norm(ck=(A*xk-b)), 
    # outs:= Struct with solver outpus
    #   outs.ex := Flag is converged to tolerance (1 conv, 0 not)
    #   outs.ctime:= Compute time used by solver
    #   out.niter:= Number of iterations
    #   out.numA:= Number of matrix multiplies A*(VEC)
    #   [out.errs]:= [Optional] If opts.store==true, then errors stored
    #   [out.times]:= [Optional] If opts.store==true, then times stored
    #-------------------------------------------------------------------------#
    # 04/02/20, J.B., Initial implementation
    # 05/24/21, J.B., Initial implementation of Kaczmark-like process
    # 05/25/21, J.B., Further modifications
    # 05/25/21, J.B., Modifications for zero residuals
    # 04/28/22, J.B., Imlementation of the extended method
    # 05/02/22, J.B., Randomized version
    # 05/03/22, J.B., Preparation to release method, to include
    #                   minimum residual and choice for permutation
    # 05/08/22, J.B., Preparation for release
    # 05/11/22, J.B., Including safeguards for the randomized method
    #                   i.e., ensuring terminating early with "zero" step
    # 05/22/22, J.B., Experiments with "zero" steps,
    #                   Python implementation 
    
    # Initializing storage    
    tstart  = time.time();
    
    (m,n) = A.shape;
    
    # Initialize storage and random permutation 
    
    #idxs    = np.random.permutation(m);
    idxs = np.arange(m);
    
    P       = np.zeros([m,n]);
    PA      = np.zeros([m,n]);
    d       = np.zeros(n);
                
    # Row index to select
    k       = 0;
    it      = 0;
    
    numA    = 0;
    ex      = 0;
    
    if store == 1:      
        
        errs  = [];
        times = [];
    
    xk      = x;
    rk      = b[idxs] - np.sum(A[idxs,:]*xk,1);
    numA    = numA + 1;
    
    nrk     = np.sqrt(np.sum((rk*rk)));
    
    # Make sure the next residual element is nonzero
    if np.abs(rk[idxs[k]]) < tol:
    
        idn = idxs[k];
    
        # Find the first instance of nonzero residual
        kn = np.nonzero(np.abs(rk[idxs])>=tol)[0];
    
        # Swap indicies
        idxs[k] = idxs[kn[0]];
        idxs[kn[0]] = idn;
        
    yk      = A[idxs[k],:];
        
    # Store minimum solution estimate
    xkmin = xk;
    rkmin = nrk;
    
    if hasMeasure == 1:
        rkm = np.sum(measure*rk,1);
        nrk = np.sqrt(np.sum((rkm*rkm)));
    else:
        nrk = np.sqrt(np.sum((rk*rk)));
    
    if store == 1:
        errs.append(nrk);
        times.append(time.time()-tstart);
    
    
    if nrk <= tol:                
        ex      = 1;
        tend    = time.time()-tstart;
        
        
    deltaki     = 1/np.sum((yk*yk));
    pk          = (deltaki*(rk[idxs[k]]))*yk;
    
    P[:,k]      = pk;
    PA[idxs,k]  = np.sum(A[idxs,:]*pk,1);
    
    d[k]        = np.sqrt(np.sum(pk*pk));
    numA        = numA + 1;
    
    if printF == 1:
    
        (m,n) = A.shape;
        
        print('');
        print('********* Linear Systems Solver ********************  ');
        print('*         Alg: PLSS KZ (Python 3.9)                   ');
        print('*         Size: m = %i, n = %i                        '% (m,n));
        print('*         tol = %1.8f                                 '% tol);
        print('*         Maxit = %i                                  '% maxiter);        
        print('****************************************************  ');        
        print('');
        
        print('Iter \t norm(Res) \t norm(pk)    ');
        print('%i \t %.4e \t %.4e   ' % (k,nrk,np.sqrt(np.sum((pk*pk)))));
        
   
    it          = it + 1;
    k           = k + 1;
    xk          = xk + pk;
    
    while np.logical_and(np.logical_and(tol < nrk,it < maxiter),ex == 0):
            
        rk[idxs] = rk[idxs] - PA[idxs,k-1];     
        
        # Make sure the next residual element is nonzero
        if np.abs(rk[idxs[k]]) < tol:
        
            idn = idxs[k];
        
            # Find the first instance of nonzero residual
            kn = np.nonzero(np.abs(rk[idxs])>=tol)[0];
        
            if np.size(kn)>0:
                # Swap indicies
                idxs[k] = idxs[kn[0]];
                idxs[kn[0]] = idn;
            else:
               ex = 1;
               tend = time.time()-tstart;
        
        yk = A[idxs[k],:];
        
        nrk     = np.sqrt(np.sum((rk*rk)));
        nyk     = np.sqrt(np.sum((yk*yk)));
                    
        rhok    = rk[idxs[k]];    
        
        
        uk1     = (PA[idxs[k],0:k])/d[0:k];
        
        nuk1    = np.sqrt(np.sum(uk1*uk1));
        
        beta    = rhok/((nyk-nuk1)*(nyk+nuk1));
        
        # Safeguard for division by zero
        if np.abs(nyk-nuk1) < 5*np.spacing(1):
            
            beta = 0.0;
            
        # Step computation        
        pk = beta*( yk - np.sum(P[:,0:k]*(uk1/d[0:k]),1) );
        
                
        # Iteration information
        if hasMeasure == 1:
            rkm = np.sum(measure*rk,1);
            nrk = np.sqrt(np.sum((rkm*rkm)));
        else:
            nrk = np.sqrt(np.sum((rk*rk))); 
            
        # Safeguard (to set pk=0 when a solution is found)
        p2      = np.sum((pk*pk));
        nrp     = np.sqrt(p2);
        if np.logical_and(nrk <= tol,np.isnan(nrp)):
            pk[:] = 0;
            nrp = 0;
                    
        if np.logical_and(printF == 1,np.logical_or(k < 11,np.mod(k,20)==0)):
            print('%i \t %.4e \t %.4e  ' % (it,nrk,nrp));         
        

        if store == 1:
            errs.append(nrk);
            times.append(time.time()-tstart);
        
        
        if nrk <= tol:        
            ex          = 1; 
            tend        = time.time()-tstart;            
            if np.logical_and(np.logical_and(printF == 1,k > 11),np.mod(k,20)!=0):
                print('%i \t %.4e \t %.4e  ' % (it,nrk,np.sqrt(np.sum((pk*pk)))));
                               
        # Store minimum iterate
        if rkmin >= nrk:            
            xkmin = xk;
            rkmin = nrk;
                    
            
        P[:,k]  = pk;
        npk     = np.sqrt(np.sum(pk*pk));
        
        if np.abs(beta) == 0.0:
            
            npk = 1.0;
        
        PA[idxs,k]  = np.sum(A[idxs,:]*pk,1);
        d[k]        = npk;
        numA        = numA + 1;
           
        # Updates or restart        
        if k == (n-1):
            k = 0;
        else:
            k = k + 1;
        
        # Prepare for next iteration
        xk  = xk + pk;
        
        it  = it + 1;
    
    if ex == 0:
        tend = time.time()-tstart;    
    
    rk = b[idxs] - np.sum(A[idxs,:]*xk,1);
    nck = np.sqrt(np.sum((rk*rk)));
    
    numA = numA+1;
    
    if nck < tol:
        ex = 1;

    
    if rkmin < nck:
       
        nrk = rkmin;
        xk = xkmin;
    
    outs = {'ex': ex,'ctime': tend, 'niter': it, 'numA': numA};
        
    if store == 1:        
        outs['errs']   = errs;
        outs['times']  = times;
    
    if printF == 1:
        
        print('');
        print('********* Summary **********************************  ');
        print('*         Conv: %i                                    ' % ex);
        print('*         Time (s): %1.2f                             ' % tend);    
        print('*         norm(Res) = %1.4f                           ' % nck);
        print('*         Iter = %i                                   ' % it);
        print('*         Num. mult. = %i                             ' % numA);
        print('****************************************************  ');
        print('');
                
    return  xk,nck,outs
