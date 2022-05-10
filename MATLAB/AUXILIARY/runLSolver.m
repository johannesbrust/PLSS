function [ outsS ] = runLSolver( A,b,x0,opts,whichSolv )
%RUNLSOLVER Runs a solver on the linear system Ax=b form
% starting point x0.
%
% INPUTS:
% A,b,x0 : Linear system inputs
% whichSolv : 'Switch' to indicate which linear solver to use
%
% OUTPUTS:
% outs.res := Exit residual
% outs.niter:= Number of iterations
% outs.time := Solver time
% outs.Info1 := Additional information
%
%--------------------------------------------------------------------------
% 08/06/20, J.B., Initial version
% 06/12/20, J.B., Scaled factorized method & Scaled factorized line search
% method
% 12/09/20, J.B., Adding "scaled recursive" PLSSR
% 05/14/21, J.B., Adding "weighted" recursive PLSS(r) solver
% 05/19/21, J.B., Adding random normal solver
% 12/06/21, J.B., Storing residual to solution (if desired) and 
%                   relative residuals
% 12/08/21, J.B., Addition of sparse random normal solver
% 01/11/22, J.B., Addition of randomized solvers
% 05/07/22, J.B., Preparation of release

optsP   = opts.optsP;
M1      = opts.M1;
M2      = opts.M2;
tolR    = opts.tolRel;
%itRes   = opts.itRes;
maxit   = opts.maxit;

ts = tic;

switch whichSolv    
    case 18 % Keep
        [xk,flag,relres,its,~] = lsqr(A,b,tolR,maxit,M1,M2,x0);       
    case 24 % Keep
        [xk,res,outs] = plss_R(x0,A,b,optsP);    
    case 28 % Keep
        [xk,res,outs] = plss_RW1(x0,A,b,optsP);
    case 29 % Keep
        [xk,res,outs] = plss_RW2(x0,A,b,optsP); % implicit ck
    case 30 % Keep
        [xk,res,outs] = plss_RANDN(x0,A,b,optsP); % random normal     
    otherwise
        [xk,res,outs] = plss_RW2(x0,A,b,optsP); % implicit ck
end

te = toc(ts);

outsS.res    = norm(A*xk-b);

outsS.time   = te;
if isfield(opts,'storeSolRes')
    if opts.storeSolRes==1
        
       outsS.solRes=norm(opts.x-xk);
        
    end
end
if isfield(opts,'storeRelRes')
    if opts.storeRelRes==1
        outsS.relRes=outsS.res/max(norm(b),1);    
    end
end

if whichSolv == 18
    outsS.niter = its;
    outsS.info1 = flag;
    outsS.info2 = relres;
else
    outsS.time  = outs.ctime;
    outsS.niter = outs.niter;
    outsS.info1 = outs;
    outsS.info2 = res;
end
