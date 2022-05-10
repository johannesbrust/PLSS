function [ xk, nck, outs ] = plss_RW2( x, A, b, opts )
%PLSS_RW2: Projected linear systems solver Algorithm 1 (Recursive) with
% inverse weighting matrix (W), (modified weighting). This implementation
% (i.e., version 2) updates residuals implicitly and explicitly forms A*pk.
%
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
% pk = -WA Sk (Sk' AWA Sk)^{-1} (Sk' Sk) e_k,                           (1)
% 
% where e_k denotes the kth column of the identity matrix. 
% When (1) is simplified, it can be expressed by the recursive formula
%
% pk = beta1*p_{k-1} + beta_2 W*yk,                                     (2)
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
%   opts.useW  := Flag to use/or not the inverse weighting
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
% 05/14/21, J.B., Weighting matrix
% 05/15/21, J.B., Modified weighting matrix
% 05/17/21, J.B., Implementation of implicit ck
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
if isfield(opts,'useW') 
    useW = opts.useW;
else
    useW = 1; 
end

% Initializing storage
tstart  = tic;

[m,n] = size(A);

% Initialize scaling matrix (and inverse)
Wh = ones(n,1);
Whi = ones(n,1);

if useW == 1
    if m ==n
        Wh(1:n) = sqrt(sum(A.*A)+1); % 0.5
        Wh(1:n) = min(Wh(1:n),1e5);
        Whi(1:n) = 1./Wh(1:n);
    else
        Wh(1:n) = 1*sqrt(sum(A.*A));    
        Whi(1:n) = 1./Wh(1:n);

        % Adjust for zero columns
        Whi(Whi==inf) = 0.0;
    end
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

yk      = A'*(ck./nck);
numA    = numA + 1;
k       = 0;

% Store minimum solution estimate
xkmin = xk;
nckmin = nck;

if hasMeasure == true; nck = norm(measure*ck); end

if nck <= tol
            
    ex      = 1;
    tend    = toc(tstart);
    
end

rhok    = nck;
zk      = Whi.*yk;
deltaki = 1/(zk'*zk);
pk      = -(deltaki*rhok).*(Whi.*zk);

if print == 1

    fprintf('\n');
    fprintf('********* Linear Systems Solver ******************** \n');
    fprintf('*         Alg: PLSS W                                \n');
    fprintf('*         Vers.: 2 (Update Res.)                     \n');
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
        
    ck      = ck + A*pk;          numA = numA + 1;    
        
    nck     = norm(ck);
    yk      = A'*(ck/nck);        numA = numA + 1;
    
    zk      = Whi.*yk;
        
    rhok    = nck;
    
    % Modifications for weighting
    Wp      = Wh.*pk;
    p2      = Wp'*Wp;    
    
    np      = sqrt(p2);    
        
    py      = pk'*yk;
    
    yy      = zk'*zk;
    
    ny      = sqrt(yy);
    
    denom   = (np*ny-py)*(np*ny+py);
    
    beta1   = (rhok * py)/denom;
    beta2   = -(rhok * p2)/denom;
    
    % Step computation
    pk      = beta1*pk + beta2*((Whi(1:n).*zk)); % .^2
        
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












