%% Demo for RandomLinearLab. Tests the overdetermined methods on a dense artificial problem
% Robert Gower and Peter Richtarik -- 18 June 2015

%------------------------- Modifications ---------------------------------%
% 04/02/20, J.B., Initial modifications for comparisons with projection 
% linear systems solvers (PLSS), J.Brust & T.Migot
%
%-------------------------------------------------------------------------%
%% Main for Probing nonlinear singular systems
clear;
close all;
%% Select Problem Parameters
m = 100;  % number of rows, 200
n =100;    %select dimension, 20
%% Sparse nonsingular
Prob = uniform_Prob(m,n);
Prob.sol = Prob.A\Prob.b;
%% Select method Parameters
options =[];
options.ep=1*10^(-4); %norm(Prob.sol)^2/10;
OUTPUTS ={};
options.max_time = 5;  % 300s
options.max_iterations = 10^6;
options.metric = eye(n); %Prob.A'*Prob.A;  % this forces all methods to use the same error measurement
%% Gaussian Least squares
[x, output_gaussls] = solve_system(Prob,@iter_Gaussls, @boot_Gaussls,options );
OUTPUTS = [ OUTPUTS ; output_gaussls];
%% LS Coordinate Descent
[x, output_lscd] = solve_system(Prob,@iter_LSCD, @boot_LSCD,options );
OUTPUTS = [ OUTPUTS ; output_lscd];
%% Gaussian Iterative Projection IGP
[x, output_gaussigp] = solve_system(Prob,@iter_Gaussigp, @boot_Gaussigp,options );
OUTPUTS = [ OUTPUTS ; output_gaussigp];
%% Randomized Karczmarz
[x, output_kz] = solve_system(Prob,@iter_Karczmarz, @boot_Karczmarz,options );
OUTPUTS = [ OUTPUTS ; output_kz];
%% plotting
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)




