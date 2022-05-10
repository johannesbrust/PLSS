%% main for testing methods for positive definite including
%  Coordinate Descent for Positive Definite  (with a block version)
%  Gaussian for Positive Definite  (with a block version)
% Robert Gower and Peter Richtarik -- 18 June 2015

%------------------------- Modifications ---------------------------------%
% 04/02/20, J.B., Initial modifications for comparisons with projection 
% linear systems solvers (PLSS), J.Brust & T.Migot
%
%-------------------------------------------------------------------------%

%% Select Dimensions
clear all
n =100;    %select dimension
% clear Prob;
%% Badly conditioned positive definites  % uncomment to test a dense positive definite matrix
% Prob = Gauss_pd_Prob(n);
%% Hilbert  % uncomment to test Hilbert matrix
% Prob.A = hilb(n);
% Prob.title =['Hilbert-' num2str(n)];
%% Sparse positive definite  % uncomment to test a sparse matrix
% density = 1/log(n^2); rc = 1/(n);  kind =1;
% Prob = sparse_symmetric_matlab(n,density,rc,kind)
%% Harwell-Boeing format
Prob.title='gr_30_30.rsa';
[ Prob.A, Prob.b] = hb_to_msm(Prob.title); 
%% Select method Parameters
options =[];
options.max_time = 5;  
options.ep=1*10^(-4);
options.max_iterations = 10^6;  % maximum of 10^6 iterations
options.metric = eye(size(Prob.A,2)); % A (default)
fix_prob_parameters; % In case only lower triangular part only, or no solution 
OUTPUTS ={};
%% Gaussian for positive definite
[x, output_gausspd] = solve_system(Prob,@iter_Gausspd, @boot_Gausspd,options );
OUTPUTS = [ OUTPUTS ; output_gausspd];
%% Randomized Coordinate Descent for positive definite
[x, output_cd] = solve_system(Prob,@iter_CD, @boot_CD,options );
OUTPUTS = [ OUTPUTS ; output_cd];
%% Block Coordinate Descent for positive definite
[x, output_Bcd] = solve_system(Prob,@iter_BCD, @boot_BCD,options );
OUTPUTS = [ OUTPUTS ; output_Bcd];
%% %% BLOCK METHODS
%% Block Gaussian positive definite
[x, output_gausspd] = solve_system(Prob,@iter_GaussBpd, @boot_GaussBpd,options );
OUTPUTS = [ OUTPUTS ; output_gausspd];
%% plotting
close all
%figure('visible','off')        
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)


