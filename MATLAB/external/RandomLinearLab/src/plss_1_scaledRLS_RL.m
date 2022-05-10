function [ xk, nck, outs ] = plss_1_scaledRLS( x, A, b, opts )
%PLSS_1_SCALEDRLS: Projected linear systems solver Algorithm 1
% Method for iteratively solving linear systems
%
% Ax = b,
%
% where A in (mxn), x in (nx1), b in (mx1) and the residuals are scaled.
%
% This algorithm computes search directions pk defined by the projections
% defined by the previous residuals in ck = Axk -b.
%
% Let Ck = [c0,...,ck] = Ck Ek Ek^{-1} = \tilde{Ck}Ek^{-1}, where
% Ek = inv(diag(norm(c0),...,norm(ck))). The step is then
%
% pk = - norm(ck)*(A \tilde{Ck} (\tilde{Ck}' AA' \tilde{Ck})^{-1} e_k),
%
% where e_k denotes the kth column of the identity matrix. Algorithm 1
% is for small/medium problems, because it stores and updates the
% inverse matrix Mk = (\tilde{Ck}' AA' \tilde{Ck})^{-1} = Rk Dk Rk'.
% This implementation includes a 'line search' parameter and an
% factorization of the middle matrix.
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
%   opts.projResTol := Tolerance for ``re-projecting"/orthgonalizing
%                       residulas.
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
% 04/06/20, J.B., Implementation of scaling version
%               , Display norm(I-Ck'*Ck,'fro') in print
% 04/09/20, J.B., Re-projecting residuals if norm(Ck(1:k-1)'ck) >
% projResTol.
% 06/08/20, J.B., Modification to number of 're-projections'
% 06/09/20, J.B., Line search
% 06/10/20, J.B., Print line search parameter
% 06/12/20, J.B., Rk Dk Rk' factorization

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
if isfield(opts,'projRes')
    projRes = opts.projRes;
else
    projRes = false;
end
if isfield(opts,'alpTol')
    alpTol = opts.alpTol
else
    alpTol = 1e-13;
end
% if isfield(opts,'projResTol')
%     projResTol = opts.projResTol;
% else
%     projResTol = 1e3;
% end

% Initializing storage
tstart  = tic;

[m,n]   = size(A);
mem     = maxiter; % Memory

Y       = zeros(n,mem);%zeros(m,mem);%Tangi
%M       = zeros(maxiter,mem);
Ck      = zeros(m,maxiter);
Imax    = eye(maxiter);
Rk      = -Imax;
Dk      = zeros(maxiter,1);
vk      = zeros(m,1);

numA    = 0;
numReproj = 0;

if store == true

    errs  = [];
    times = [];

end

xk      = x;
Ax      = A*xk;
ck      = Ax-b;
%ck      = A*xk-b;
rhok    = norm(ck);
%rhokp   = rhok;
numA    = numA + 1;
cks     = ck./rhok;
yk      = A'*cks; % Scaling method
numA    = numA + 1;
k       = 0;

if hasMeasure == true; nck = norm(measure*ck); else nck = rhok; end

% Initializing minimum residual iterate
xkmin    = xk;
nckmin  = nck;
kmin    = k;

if nck <= tol

    tend        = toc(tstart);
    outs.ex     = 1;
    outs.ctime  = tend;
    outs.niter  = k;
    outs.numA   = numA;
    outs.numReproj = numReproj;
    outs.kmin = kmin;

    if store == true
        errs        = [errs;nck];
        times       = [times;tend];
        outs.errs   = errs;
        outs.times  = times;
    end

    return;

end

%rhok        = ck'*ck;
deltaki     = 1/(yk'*yk);
pk          = -(deltaki*rhok).*yk;

Ck(:,k+1)   = cks;

%alphak      = rhok^2/norm(A*pk)^2;
Ap          = A*pk;
numA        = numA + 1;
alphak      = rhok^2/(Ap'*Ap);

if print == 1

    fprintf('----------- Running Algorithm 1 ----------- \n'); % 43 chars
    fprintf('Iter \t norm(ck)      \t norm(pk) \t err(CkTCk-I) \t alphak  \n');

    errorth = norm(Ck(:,1)'*Ck(:,1)-Imax(1,1),'fro');
    fprintf('%i \t %.4e \t %.4e \t %.4e \t %.4e \n',k,nck,norm(pk),errorth,alphak);

end

k           = k + 1;
Y(:,k)      = yk;
%M(1:k,1:k)  = deltaki;%1/(yk'*yk);
Dk(1)       = deltaki;

Ax          = Ax + alphak.*Ap;
xk          = xk + alphak.*pk;
%xk          = xk + pk;

% TMP timings
%cTimes = zeros(3,1);

while (tol < nck) && (k < maxiter) %&& (toc(tstart) < maxtime)

    % A applies and 'reprojection'
    %tA = tic;

    %ck      = A*xk - b;             numA = numA + 1;
    ck      = Ax - b;

    rhok    = norm(ck);
    cks     = ck./rhok;

    % ``Re-projection"
    ck1 = Ck(:,1:k)'*ck;


     %ortherr = norm(Ckm1ck);

     if projRes == true
        Ckm1ck  = ck1/rhok;
%    if alphak < projResTol
        cks = cks - Ck(:,1:k)*((Ck(:,1:k)'*Ck(:,1:k))\Ckm1ck);
        %cks = cks - (Ck(:,1:k)*((Ck(:,1:k)'*Ck(:,1:k))\Ckm1ck)/rhok);
        %cks = cks - (1-alphak)*(rhokp/rhok).*Ck(:,k);
        %cks = cks - Ck(:,1:k)*Ckm1ck;
        %cks = cks
        numReproj = numReproj + 1;
    end

    yk      = A'*cks;        numA = numA + 1;

 %   cTimes(1) = cTimes(1) + toc(tA);
    %rhok    = ck'*ck;

    % Step
   % tS = tic;

    uk1     = Y(:,1:k)'*yk;
    %uk      = M(1:k,1:k)*uk1;
    uk      = Rk(1:k,1:k)*(Dk(1:k).*(Rk(1:k,1:k)'*uk1));

    deltaki = 1/(yk'*yk - uk'*uk1);

    Rk(1:k,k+1) = uk;
    Dk(k+1)     = deltaki;

    %vk          = Rk(1:k,1:k)*(Dk(1:k).*(Rk(1:k,1:k)'*Ckm1ck));
    %gammak      = uk'*Ckm1ck - rhok^2;

    % Step computation
    %pk      =-(rhok*deltaki).*(yk-Y(:,1:k)*uk(1:k));

 %   cTimes(2) = cTimes(2) + toc(tS);

    % Update middle matrix
    %tM = tic;

    %M(1:k,1:k) = M(1:k,1:k) + (deltaki.*uk(1:k))*uk(1:k)';
    %M(1:k,k+1) = -deltaki.*uk(1:k);
    %M(k+1,1:k) = M(1:k,k+1)';%-deltaki.*uk(1:k)';
    %M(k+1,k+1) = deltaki;
    %M(1:(k+1),1:(k+1)) = deltaki.*M(1:(k+1),1:(k+1));

 %   cTimes(3) = cTimes(3) + toc(tM);

    Y(:,k+1)    = yk;
    Ck(:,k+1)   = cks;
    vk(1:k,1)   = ck1; %[ck1;cks'*ck]; % Ck(:,1:k+1)'*ck; %
    vk(k+1,1)   = cks'*ck;
    % Step computation
    %pk      =-Y(:,1:k+1)*M(1:k+1,1:k+1)*(Ck(:,1:k+1)'*ck);
    pk      = -Y(:,1:k+1)*(Rk(1:k+1,1:k+1)*(Dk(1:k+1).*(Rk(1:k+1,1:k+1)'*vk(1:k+1,1))));

    % Iteration information
    if hasMeasure == true; nck = norm(measure*ck); else nck = norm(ck); end

    if nck < nckmin
        xkmin = xk;
        nckmin = nck;
        kmin = k;
    end

    Ap          = A*pk;
    numA        = numA + 1;
    alphak      = rhok^2/(Ap'*Ap);
    %alphak      = rhok^2/norm(A*pk)^2;
    if print == 1

        errorth = norm(Ck(:,1:(k+1))'*Ck(:,1:(k+1))-Imax(1:(k+1),1:(k+1)),'fro');
        fprintf('%i \t %.4e \t %.4e \t %.4e \t %.4e \n',k,nck,norm(pk),errorth,alphak);

        %fprintf('%i \t %.4e \t %.4e \n',k,nck,norm(pk));
    end
    if store == true
        errs        = [errs;nck];
        times       = [times;toc(tstart)];
    end
    
    if nck <= tol

        tend        = toc(tstart);
        outs.ex     = 1;
        outs.ctime  = tend;
        outs.niter  = k;
        outs.numA   = numA;
        outs.numReproj = numReproj;

        if store == true
            outs.errs   = errs;
            outs.times  = times;
        end

        return;
    end
    
    if alphak < alpTol
        break;
    end

    % Prepare for next iteration
%     Ap          = A*pk;
%     numA        = numA + 1;
%     alphak      = rhok^2/(Ap'*Ap);

    Ax          = Ax + alphak.*Ap;
    xk          = xk + alphak.*pk;

    %xk          = xk + alphak.*pk;
    %xk  = xk + pk;
    %rhokp = rhok;
    k   = k + 1;
end

ck   = A*xk - b;
numA = numA + 1;
if hasMeasure == true; nck = norm(measure*ck); else nck = norm(ck); end

if nckmin < nck
   xk  = xkmin;
   nck = nckmin;
   outs.kmin = kmin;
else
   outs.kmin = k;
end

tend        = toc(tstart);
if nckmin < tol; outs.ex = 1;else outs.ex =0; end
outs.ctime  = tend;
outs.niter  = k;
outs.numA   = numA;
outs.numReproj = numReproj;

if store == true
    outs.errs   = errs;
    outs.times  = times;
end

end
