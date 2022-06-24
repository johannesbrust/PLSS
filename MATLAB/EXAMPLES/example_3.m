%------------------------ example_3 --------------------------------------%
% 
% Example of applying the generalized Kaczmarz process to an
% underdetermined system
% The matrix is build-in in Matlab
%
% m = 20; n = 1000; N = [m,n];
% A = gallery('kahan',N);
%
%-------------------------------------------------------------------------%
% 06/24/22, J.B., Preparation for release

clc;
clear;

addpath('../ALGS');

% Setup problem
m       = 200; n = 1000; N = [m,n];
A       = gallery('kahan',N);

b       = full(sum(A,2));
nb      = norm(b);
x0      = zeros(size(A,2),1);

% Solver options
opts.print = 1;
opts.store = 1;
opts.maxiter = size(A,1)+10;

% Call solver
[xk,nck,outs] = plss_KZ_EX(x0,A,b,opts);