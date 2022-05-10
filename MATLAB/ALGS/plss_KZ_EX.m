function [ xk, nck, outs ] = plss_KZ_EX( x, A, b, opts )
%PLSS_KZ_EX: Projected linear systems solver using a Kaczmark like process. 
%
% Method for iteratively solving linear systems
% 
% Ax = b,
%
% where A in (mxn), x in (nx1), b in (mx1).
%
% This algorithm computes search directions dk defined by the projections
% defined by sketches that store previous residuals in -ck = rk = b - A*xk.
% 
% Let Sk = [e1,...,ek1], then the step is
% 
% pk = -A Sk (Sk' AA Sk)^{-1} (Sk' ck),
% 
% where e_k denotes the kth column of the identity matrix. This method
% stores the updates Pk=[p0,...,pkm1] and the products PkA =
% [A*p0,...,A*pkm1]. 
% Options for permutation of rows is possible:
%       perm =  0 : No reordering
%               1 : Increasing row norms
%               2 : Decreasing residuals
%               3 : dmperm "Dulmage-Mendelsohn permutation"
%               4 : randomized (default)
% Note, on restart the iteration starts with decreasing residuals
%
% INPUTS:
% x := Initial guess (nx1)
% A := Linear system (mxn)
% b := Right hand side (mx1)
% opts := Struct with solver options
%   opts.tol    := Residual tolerance 
%   opts.maxiter  := Maximum iterations
%   opts.maxtime  := Maximum time
%   opts.measure:= Different error measure (for comparisons). 
%   opts.print  := Flag for printing
%   opts.store  := Flag for storing information (for comparisons)
%   opts.rowTol  := Row norm tolerance
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
% 05/24/21, J.B., Initial implementation of Kaczmark-like process
% 05/25/21, J.B., Further modifications
% 05/25/21, J.B., Modifications for zero residuals
% 04/28/22, J.B., Imlementation of the extended method
% 05/02/22, J.B., Randomized version
% 05/03/22, J.B., Preparation to release method, to include
%                   minimum residual and choice for permutation
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
if isfield(opts,'rowTol') 
    rowTol = opts.rowTol;
else
    rowTol = 1e-12; 
end
if isfield(opts,'perm') 
    perm = opts.perm;
else
    perm = 4; 
end

% Initializing storage
tstart  = tic;

[m,n]   = size(A);
%mem     = maxiter; % Memory

% Sparse storage
P       = sparse(m,n);
PA      = sparse(m,n);
d       = zeros(n,1);

numA    = 0;
ex      = 0;

if store == true
   
    errs  = [];
    times = [];
    
end

% Row index to select
k       = 0;
it      = 0;

% Initializing permutation strategy
if perm == 0
    
    % Original rows selected 
    idxs = 1:m;
elseif perm == 1
    
    % Inspect matrix for zero rows,
    % by computing row norms "rn". Selection of rows is
    % in increasing norm magnitutes
    AA = A.*A;
    rn = sqrt(sum(AA,2));

    [srn,idxs1] = sort(rn);
    
    % Search for index such that row norm is nonzero
    rnidx = 1;
    for i = 1:m
        if srn(i) > rowTol
            break;
        else
            rnidx = rnidx + 1;
        end
    end
    idxs = idxs1(rnidx:end);
    
elseif perm == 2
    
    % Use residuals in decreasing order
    [~,idxs] = sort(abs(b),'descend');
    
elseif perm == 3
    
    % Dmperm function call
    idxs = dmperm(A);
    
elseif perm == 4
    
    % Random permutation
    idxs = randperm(m);
    
end

xk      = x;
rk      = b(idxs)-A(idxs,:)*xk;
numA    = numA + 1;

yk      = A(idxs(k+1),:)'; 

nck     = norm(rk);

% Store minimum iteration
xkmin   = xk;
nckmin  = nck;

if hasMeasure == true; nck = norm(measure*ck); end % else nck = norm(ck); end

if store == true
    errs        = [errs;nck];
    times       = [times;toc(tstart)];    
end

if nck <= tol
            
    ex      = 1;
    tend    = toc(tstart);
    
end

deltaki     = 1/(yk'*yk);
pk          = (deltaki*rk((1:k+1))).*yk;

P(:,k+1)    = pk;
PA(:,k+1)   = A(idxs,:)*pk;
d(k+1)      = sqrt(pk'*pk);
numA        = numA + 1;

if print == 1

    fprintf('\n');
    fprintf('********* Linear Systems Solver ******************** \n');
    fprintf('*         Alg: PLSS KZ                               \n');    
    fprintf('*         Size: m = %i, n = %i                       \n',m,n);
    fprintf('*         tol = %1.8f                                \n',tol);
    fprintf('*         Maxit = %i                                 \n',maxiter);
    fprintf('**************************************************** \n');
    fprintf('\n');
        
    fprintf('Iter \t norm(Res) \t norm(pk)   \n');
    fprintf('%i \t %.4e \t %.4e  \n',k,nck,norm(pk));
    
end

it  = it + 1;
k   = k + 1;
kk  = 1;
xk  = xk + pk;

% Main loop
while (tol < nck) && (it < maxiter) && (ex == 0)
    
    rk      = rk - PA(:,k); 
        
    yk = A(idxs(kk+1),:)';
    
    nck = norm(rk);
    nyk = norm(yk);
    
    rhok = rk((kk+1)); 
        
    uk1     = ((PA(kk+1,1:k))')./d(1:k);
    nuk1    = norm(uk1);
    beta    = rhok/((nyk-nuk1)*(nyk+nuk1));
    
    % Step computation
    pk = beta.*(yk - P(:,1:k)*(uk1./d(1:k)));
               
   % Iteration information
    if hasMeasure == true; nck = norm(measure*rk); end    % else nck = norm(ck); 
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
        
    % Updates or restart
    if mod(it,(n-1))~=0
        npk = norm(pk);
        if npk > 1e-14
            P(:,k+1)    = pk; %#ok<SPRIX>
            PA(:,k+1)   = A(idxs,:)*pk; %#ok<SPRIX>
            d(k+1)      = npk;%sqrt(pk'*pk);
            numA        = numA + 1;
            k           = min(k + 1,n-1);
        end
    else
        % Restart
        % Including new sort of residuals
        k   = 0;
        kk  = 0;
        rk  = b(idxs)-A(idxs,:)*xk;
        
        [~,idxs]    = sort(abs(rk),'descend');
        
        yk          = A(idxs(k+1),:)'; %'*ck;
        numA        = numA + 1;

        %rhok        = ck'*ck;
        deltaki     = 1/(yk'*yk);
        pk          = (deltaki*rk((1:k+1))).*yk;

        P(:,k+1)    = pk; %#ok<SPRIX>
        PA(:,k+1)   = A(idxs,:)*pk; %#ok<SPRIX>
        d(k+1)      = sqrt(pk'*pk);
        numA        = numA + 1;
        k           = min(k + 1,n-1);
        
    end
        
    it      = it +1;    
    kk      = min(kk + 1,n-1);
    
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















