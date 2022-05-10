% main for testing optimized probabilitites of the Kaczmarz method
% Requires cvx  see http://cvxr.com/cvx/
% 2015 - Robert. M Gower and Peter Richtárik
cvx_setup
%% main for positive definite systems
n =50;    %select dimension
clear Prob;
%% Badly conditioned positive definites
Prob = Gauss_pd_Prob(n);
%% Parameters
options =[];
options.max_time = 10;  % 10 secs
options.ep=1*10^(-4); %norm(Prob.sol)^2/10;
options.max_iterations = sum(size(Prob.A))^3;
fix_prob_parameters;
OUTPUTS ={};
%% Calculate optimal CD prob
 Prob.title = [ Prob.title '-opt'];
display( Prob.title);
n = length(Prob.b);
[p, Value_OPT] =  opt_rates(Prob.A,Prob.A,eye(n,n));
%% Randomized Coordinate Descent for positive definite
[x, output_cd] = solve_system(Prob,@iter_CD, @boot_CD,options );
OUTPUTS = [ OUTPUTS ; output_cd];
%% Randomized Coordinate Descent for positive definite with user prob
options.probs = p;
[x, output_cd] = solve_system(Prob,@iter_CD, @boot_CD,options );
OUTPUTS = [ OUTPUTS ; output_cd];
    %% plotting
close all     
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)



