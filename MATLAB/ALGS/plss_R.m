function [ xk, nck, outs ] = plss_R( x, A, b, opts )
%PLSS_R: Projected linear systems solver Algorithm 1 (Recursive)
% Method for iteratively solving linear systems
% 
% Ax = b,
%
% where A in (mxn), x in (nx1), b in (mx1).
%
% This algorithm computes search directions pk defined by the projections
% defined by sketches that store previous residuals in -ck = rk = b - A*xk.
% 
% Let Sk = [c0,...,ck], then the step is
% 
% pk = -A Sk (Sk' AA Sk)^{-1} (Sk' Sk) e_k,                             (1)
% 
% where e_k denotes the kth column of the identity matrix. 
% When (1) is simplified, it can be expressed by the recursive formula
%
% pk = beta1*p_{k-1} + beta_2*yk,                                       (2)
%
% where yk = A'ck.
%
% INPUTS:
% x := Initial guess (nx1)
% A := Linear system (mxn)
% b := Right hand side (mx1)
% opts := Struct with solver options
%   opts.tol    := Residual tolerance 
%   opts.maxiter  := Maximum iterations
%   opts.measure:= Different error measure (for comparisons). 
%   opts.print  := Flag for printing
%   opts.store  := Flag for storing information (for comparisons)
%
% OUTPUTS:
% xk := Iterate at termination
% nck:= norm(ck=(A*xk-b)), 
% outs:= Struct with solver outpus
%   outs.ex := Flag is converged to tolerance (1 conv, 0 not)
%   outs.ctime:= Compute time used by solver
%   out.niter:= Number of iterations
%   out.numA:= Number of matrix multiplies A*(VEC)
%   [out.errs]:= [Optional] If opts.store==true, then errors stored
%   [out.times]:= [Optional] If opts.store==true, then times stored
%-------------------------------------------------------------------------%
% 04/02/20, J.B., Initial implementation
% 12/07/20, J.B., Recursive formula
% 05/08/22, J.B., Preparation for release
% 11/08/22, J.B., Safeguarding residual norm on exit

% Initializations
if isfield(opts,'tol') 
    tol = opts.tol;
else
    tol = 1e-6; 
end
if isfield(opts,'maxiter') 
    maxiter = opts.maxiter;
else
    maxiter = 500; 
end
% if isfield(opts,'maxtime') 
%     maxtime = opts.maxtime;
% else
%     maxtime = 5; 
% end
if isfield(opts,'measure') 
    measure = opts.measure;
    hasMeasure = true;
else    
    hasMeasure = false;
end
if isfield(opts,'print') 
    print = opts.print;
else
    print = 0; 
end
if isfield(opts,'store') 
    store = opts.store;
else
    store = 0; 
end

% Initializing storage
tstart  = tic;

numA    = 0;
ex      = 0;
if store == true
   
    errs  = [];
    times = [];
    
end

xk      = x;
ck      = A*xk-b;
numA    = numA + 1;

nck     = norm(ck);
yk      = A'*(ck/nck);

numA    = numA + 1;
k       = 0;

% Store minimum solution estimate
xkmin = xk;
nckmin = nck;

if hasMeasure == true; nck = norm(measure*ck); end % else nck = norm(ck); end

if store == true
    errs        = [errs;nck];
    times       = [times;toc(tstart)];    
end

if nck <= tol
            
    ex      = 1;
    tend    = toc(tstart);
        
    outs.ex     = ex;
    outs.ctime  = tend;
    outs.niter  = k;
    outs.numA   = numA;
    
    if store == true        
        outs.errs   = errs;
        outs.times  = times;
    end
    
    return;
    
end

rhok        = nck;
deltaki     = 1/(yk'*yk);
pk          = -(deltaki*rhok).*yk;

if print == 1

    [m,n] = size(A);
    
    fprintf('\n');
    fprintf('********* Linear Systems Solver ******************** \n');
    fprintf('*         Alg: PLSS                                  \n');    
    fprintf('*         Size: m = %i, n = %i                       \n',m,n);
    fprintf('*         tol = %1.8f                                \n',tol);
    fprintf('*         Maxit = %i                                 \n',maxiter);
    fprintf('**************************************************** \n');
    fprintf('\n');
        
    fprintf('Iter \t norm(Res) \t norm(pk)   \n');
    fprintf('%i \t %.4e \t %.4e  \n',k,nck,norm(pk));
    
end

k           = k + 1;
xk          = xk + pk;

while (tol < nck) && (k < maxiter) && (ex == 0)
    
    ck      = A*xk - b;     numA = numA + 1;    
    
    nck     = norm(ck);
    
    yk      = A'*(ck/nck);  numA = numA + 1;
        
    rhok    = nck;    
    
    p2      = pk'*pk;
    np      = sqrt(p2);    
    
    py      = pk'*yk;
    yy      = yk'*yk;
    ny      = sqrt(yy);
    
    denom   = (np*ny-py)*(np*ny+py);
    
    beta1   = (rhok * py)/denom;
    beta2   = -(rhok * p2)/denom;
    
    % Step computation
    pk      = beta1*pk + beta2*yk;
        
    % Iteration information
    if hasMeasure == true; nck = norm(measure*ck); end    % else nck = norm(ck); 
    if print == 1 && ((k < 11) || mod(k,20)==0)
        fprintf('%i \t %.4e \t %.4e \n',k,nck,norm(pk)); 
    end
    
    if store == true
        errs        = [errs;nck]; %#ok<AGROW>
        times       = [times;toc(tstart)]; %#ok<AGROW>
    end
    
    if nck <= tol        
        ex          = 1; 
        tend        = toc(tstart);            
        if print == 1 && 11 < k && mod(k,20)~=0
            fprintf('%i \t %.4e \t %.4e \n',k,nck,norm(pk));
        end
    end    
    
    % Store minimum iterate
    if nckmin >= nck
        
        xkmin = xk;
        nckmin = nck;
        
    end
    
    % Prepare for next iteration
    xk  = xk + pk;
    k   = k + 1;
end

if ex == 0
    tend = toc(tstart);    
end

nck = norm(A*xk-b);

if nck < tol
    ex = 1;
end

if nckmin < nck || isnan(nck)
   
    nck = nckmin;
    xk = xkmin;
    
end

outs.ex     = ex;
outs.ctime  = tend;
outs.niter  = k;
outs.numA   = numA;
    
if store == true        
    outs.errs   = errs;
    outs.times  = times;
end

if print == 1
    
    fprintf('\n');
    fprintf('********* Summary ********************************** \n');
    fprintf('*         Conv: %i                                   \n',ex);
    fprintf('*         Time (s): %1.2f                            \n',tend);    
    fprintf('*         norm(Res) = %1.4f                          \n',nck);
    fprintf('*         Iter = %i                                  \n',k);
    fprintf('*         Num. mult. = %i                            \n',numA);
    fprintf('**************************************************** \n');
    fprintf('\n');
            
end
    














