%------------------------ example_1 --------------------------------------%
% 
% Example of applying PLSS (Algorithm 1, W=I)
% to a sparse Poisson linear systems 
% The matrix is build-in in Matlab
%
% n1 = 20; n = n1^2;
% A = gallery('poisson',n1);
%
% The solution is x = ones(n,1);
%-------------------------------------------------------------------------%
% 05/08/22, J.B., Preparation for release

clc;
clear;

addpath('../ALGS');

% Setup problem
n1      = 20;
A       = gallery('poisson',n1);
[m,n]   = size(A);

b       = full(sum(A,2));
nb      = norm(b);
x0      = zeros(size(A,2),1);

% Solver options
opts.print = 1;
opts.store = 1;
opts.maxiter = size(A,2);

% Call solver
[xk,nck,outs] = plss_R(x0,A,b,opts);

