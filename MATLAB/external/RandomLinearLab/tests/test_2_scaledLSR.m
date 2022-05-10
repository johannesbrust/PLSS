%------------------------ test_2_scaledLS --------------------------------%
% 
% Test using a Matlab, cgs and minres solver and the scaled
% small/medium plss (Alg. 1)
%-------------------------------------------------------------------------%
% 04/02/20, J.B., Initial version
% 04/03/20, J.B., Including cgs, minres
% 04/06/20, J.B., Scaled version
% 06/09/20, J.B., Scaled Y version, LS version
% 06/12/20, J.B., Scaled Y version with factorization
clc;
clear;

addpath('../src');

% Using the example from Matlab's documentation on gmres
load west0479;
A   = west0479;
b   = full(sum(A,2));
nb  = norm(b);
x0  = zeros(size(A,2),1); 
[m,n] = size(A);

% Solver options
opts.print = 1;
%opts.tol = 1e-6*nb;

opts.projRes = true;
opts.maxiter    = m;
%opts.print = 0;
[xksYR,ncksYR,outssYR] = plss_1_scaledRLS(x0,A,b,opts);

% cgs 'Conjugate Gradients Squared Method'
% Uses normalized residual to stop
tol     = 1e-6/nb;
maxit   = m;
M1      = []; % preconditioner 1
M2      = []; % preconditioner 2
[X,FLAG,RELRES,ITER,RESVEC] = cgs(A,b,tol,maxit,M1,M2,x0);

% gmres
tgm = tic;
[xg,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,M1,M2,x0);
tgme = toc(tgm);