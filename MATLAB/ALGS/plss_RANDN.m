function [ xk, nck, outs ] = plss_RANDN( x, A, b, opts )
%PLSS_RANDN: Projected linear systems solver with full random 
% normal matrix
% Method for iteratively solving linear systems
% 
% Ax = b,
%
% where A in (mxn), x in (nx1), b in (mx1).
%
% This algorithm computes search directions dk defined by the projections
% defined by the random matrix Sk = randn(n,r) where "r" is a paramter
%  
% pk = - A Sk (Sk' AA Sk)^{-1} (Sk' (Axk-b))
% 
% INPUTS:
% x := Initial guess (nx1)
% A := Linear system (mxn)
% b := Right hand side (mx1)
% opts := Struct with solver options
%   opts.tol    := Residual tolerance 
%   opts.maxiter  := Maximum iterations
%   opts.measure:= Different error measure (for comparisons)
%   opts.print  := Flag for printing
%   opts.store  := Flag for storing information (for comparisons)
%   opts.r  := Columns of the random sketching matrix
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
% 05/19/21, J.B., Implementation of random matrices
% 05/08/22, J.B., Preparation for release

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

[m,n]   = size(A);

if isfield(opts,'r') 
    r = opts.r;
else
    r = ceil(min(50,0.1*n));
end

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

Ck = randn(m,r);
Yk = A'*Ck;
yk = Ck'*ck;

numA = numA + r;
k       = 0;

% Store minimum solution estimate
xkmin = xk;
nckmin = nck;

if hasMeasure == true; nck = norm(measure*ck); end

if store == true
    errs        = [errs;nck];
    times       = [times;toc(tstart)];    
end

if nck <= tol
            
    ex      = 1;
    tend    = toc(tstart);
    
end

pk          = -Yk*((Yk'*Yk)\yk);

if print == 1

    fprintf('\n');
    fprintf('********* Linear Systems Solver ******************** \n');
    fprintf('*         Alg: Random Normal Sketch                  \n');    
    fprintf('*         Size: m = %i, n = %i                       \n',m,n);
    fprintf('*         tol = %1.8f                                \n',tol);
    fprintf('*         Maxit = %i                                 \n',maxiter);
    fprintf('*         r = %i                                     \n',r);
    fprintf('**************************************************** \n');
    fprintf('\n');
        
    fprintf('Iter \t norm(Res) \t norm(pk)   \n');
    fprintf('%i \t %.4e \t %.4e  \n',k,nck,norm(pk));
    
end

k           = k + 1;
xk          = xk + pk;

while (tol < nck) && (k < maxiter) && (ex == 0)
    
    ck      = A*xk - b;    numA = numA + 1;    
    
    nck     = norm(ck);
    
    Ck = randn(m,r);
    Yk = A'*Ck;
    yk = Ck'*ck;

    numA = numA + r;

    pk = -Yk*((Yk'*Yk)\yk);
        
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

if nckmin < nck
   
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














