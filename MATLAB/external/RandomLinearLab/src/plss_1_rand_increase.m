function [ xk, nck, outs ] = plss_1_rand_increase( x, A, b, opts, Prob )
%PLSSR: Projected linear systems solver Algorithm 1 randomized
% Method for iteratively solving linear systems
% 
% Ax = b,
%
% where A in (mxn), x in (nx1), b in (mx1).
%
% This algorithm computes search directions dk defined by the projections
% defined by the previous residuals in ck = Axk -b.
% 
% Let Ck = [c0,...,ck], then the step is
% 
% Yk = A Ck
% pk = - Yk (Yk' Yk)^{-1} (Ck' Ck) e_k,
% 
% where e_k denotes the kth column of the identity matrix. 
% The limited memory parameter (mem) limit the number of residuals stored.
% We store scaled residuals
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
%   opts.zeroStart := Flag if starting point is zero
%   opts.mem    := Number of residuals stored
%   opts.rst    := Restart every rst iterations.
%
% OUTPUTS:
% xk := Iterate at termination
% nck:= norm(ck=(A*xk-b)), 
% outs:= Struct with solver outpus
%   outs.ex     := Flag is converged to tolerance (1 conv, 0 not)
%   outs.ctime  := Compute time used by solver
%   out.niter   := Number of iterations
%   out.numA    := Number of matrix multiplies A*(VEC)
%   out.nRst    := Number of restart
%   [out.errs]  := [Optional] If opts.store==true, then errors stored
%   [out.times] := [Optional] If opts.store==true, then times stored
%-------------------------------------------------------------------------%
% 04/02/20, J.B., Initial implementation
% 07/15/20, T.M., randomized method selecting a random fixed number of rows of A.
% 09/06/20, T.M., randomized method selecting an increasing number of rows of A.
% 09/06/20, T.M, real tracking of the norm of the residual (not the estimated one)
                 

% Initializations
if isfield(opts,'tol') 
    tol = opts.tol;
else
    tol = 1e-6; 
end
if isfield(opts,'maxiter') 
    maxiter = opts.maxiter;
else
    maxiter = 20; 
end
if isfield(opts,'maxtime') 
    maxtime = opts.maxtime;
else
    maxtime = 5; 
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
if isfield(opts,'zeroStart')
    zS = opts.zeroStart;
else
    zS = 0; %1 if zeros(n,1) is starting point.
end
if isfield(opts,'initsampleSize')
    splsze = opts.initsampleSize; %initial size of the sampleSize.
else
    splsze = 2;
end
if isfield(opts,'rst')
    rst = opts.rst;
else
    rst = 250;
end
if isfield(opts,'mem') %Memory
     mem = opts.mem;
else
     mem = splsze; %or maxiter
end
if isfield(opts,'errorScale') 
    errorScale = opts.errorScale;
else
    errorScale = 100; 
end
if isfield(opts,'errorSqrt') 
    errorSqrt = opts.errorSqrt;
else
    errorSqrt = 0; 
end

% Initializing storage
tstart  = tic;

[m,n]   = size(A);

Ck      = zeros(m,mem);

splsze     = min(splsze, m);
max_splsze = n;

numA    = 0;
mIdx    = 1:mem;
nMem    = 0;
nRst    = 0;

if store == true
   
    errs  = [];
    times = [];
    
end

xk      = x;
if zS == 1; ck = -b; else ck = A*xk-b; numA = numA + 1; end
k       = 0;

if hasMeasure == true; nck = norm(measure*ck); else nck = norm(ck); end

if nck <= tol
        
    tend        = toc(tstart);    
    outs.ex     = 1;
    outs.ctime  = tend;
    outs.niter  = k;
    outs.numA   = numA;
    outs.nRst   = nRst;
    
    if store == true
        errork      = errorScale*calculate_error(x,Prob,opts);
        errs        = [errs;errork ];
        times       = [times;tend];
        outs.errs   = errs;
        outs.times  = times;
    end
    
    return;
    
end

rdi         = sort(randperm(m, splsze)); %randperm(m, splsze);%randperm(m, k+1);
Ardi        = A(rdi,:);
ckrdi       = ck(rdi);
nckrdi      = norm(ckrdi);
ckrdis      = ckrdi ./ nckrdi;
yk          = Ardi' * ckrdi;            numA        = numA + 1;
deltaki     = 1 / norm(yk)^2;
pk          = -deltaki * yk; %-(deltaki*rhok).*yk; %scaled residuals

if print == 1

    fprintf('----------- Running Algorithm 1 ----------- \n'); % 43 chars
    fprintf('Iter \t norm(ck)      \t norm(pk) \t norm(cki)   \n');
    fprintf('%i \t %.4e \t %.4e \t %.4e  \n',k,nck,norm(pk),norm(ckrdi)*m/splsze);
    
end

nMem               = nMem+1;
k                  = k + 1;
Ck(rdi,mIdx(nMem)) = ckrdis;
xk                 = xk + pk;


while (tol < nck) && (k < maxiter) && (toc(tstart) < maxtime)

    splsze  = min(splsze + 1,max_splsze);
    rdi     = sort(randperm(m, splsze));
    Ardi    = A(rdi,:);
    ckrdi   = Ardi * xk - b(rdi);
    nckrdi  = norm(ckrdi);
    ckrdis  = ckrdi ./ nckrdi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ck      = A*xk - b;    numA = numA + 1;    
    %yk      = A'*ck;        numA = numA + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nMem < mem
        nMem = nMem+1;
        Ck(rdi,mIdx(nMem))   = ckrdis;
    else
        ms               = mIdx(1);
        Ck(rdi,ms)       = ckrdis;
        mIdx(1:(nMem-1)) = mIdx(2:nMem);
        mIdx(nMem)       = ms;
    end

    Yk          = Ardi' * Ck(rdi,1:nMem);                  numA        = numA + 1;
    %Mk          = inv(Yk' * Yk + tol * eye(nMem));
    %uk          = Mk * (Ck(rdi,1:nMem)' * ckrdi);
    uk          = (Yk' * Yk + tol * eye(nMem)) \ (Ck(rdi,1:nMem)' * ckrdi);
    %uk          = cgs((Yk' * Yk + tol * eye(nMem)),(Ck(rdi,1:nMem)' * ckrdi),tol);
    pk          = - Yk * uk;

    % Iteration information
    oldnck = nck;
    nck    = nckrdi*m/splsze; %average estimation of the residual   
    xk     = xk + pk;

    if print == 1
        fprintf('%i \t %.4e \t %.4e \t %.4e \n',k,norm(ck),norm(pk),nck); 
    end   
    if store == true
        err         = calculate_error(xk,Prob,opts);
        errork      = errorScale*err;
        errs        = [errs;errork ];
        times       = [times;toc(tstart)];
        if(err < opts.ep)
        	fprintf('%d  | %3.4f  |  %3.4f \n',k,errork,times(k) );
        	break;
    	end
    end
    
    if nck <= tol
        
        tend        = toc(tstart);    
        outs.ex     = 1;
        outs.ctime  = tend;
        outs.niter  = k;
        outs.numA   = numA;
        outs.nRst   = nRst;

        if store == true        
            outs.errs   = errs;
            outs.times  = times;
        end

        return;
    end    
    
    % Restart
    if mod(k,rst) == 0
        nMem = 0;
        mIdx = 1:mem;
        nRst = nRst+1;
    end

    % Prepare for next iteration
    %[k, nck, opts.maxiter]
    %xk  = xk + pk;
    k   = k + 1;
end

tend        = toc(tstart);    
outs.ex     = 0;
outs.ctime  = tend;
outs.niter  = k;
outs.numA   = numA;

if store == true        
    outs.errs   = errs;
    outs.times  = times;
end

end %end of function













