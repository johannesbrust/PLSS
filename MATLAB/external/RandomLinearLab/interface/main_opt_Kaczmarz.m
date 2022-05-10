% main for testing optimized probabilitites of the Kaczmarz method
% Requires cvx  see http://cvxr.com/cvx/
% 2015 - Robert. M Gower and Peter Richtárik
%% main for positive definite systems
clear all
m = 500;  % number of rows
n =100;    %select dimension
%% Badly conditioned positive definites
Prob = uniform_Prob(m,n);
%% Calculate optimal Kaczmarz prob
display([ Prob.title '-popt-k'])
[p, Value_OPT,optp_time] =  opt_rates(Prob.A,eye(size(Prob.A,2)),eye(size(Prob.A,1)));
 display(optp_time);
%% Randomized Karczmarz
% preparing row-weighted probability dist
[x, output_kz] = solve_system(Prob,@iter_Karczmarz, @boot_Karczmarz,options );
OUTPUTS = [ OUTPUTS ; output_kz];
%% Optimized Randomized Karczmarz
% preparing row-weighted probability dist
options.probs = p;
[x, output_kz] = solve_system(Prob,@iter_Karczmarz, @boot_Karczmarz,options );
OUTPUTS = [ OUTPUTS ; output_kz];
Prob.title = [ Prob.title '-popt-k']
%% plotting
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)