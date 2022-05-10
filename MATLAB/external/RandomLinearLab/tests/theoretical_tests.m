%% Theoretical convergence tests
n=50;
m =500; clear Prob;
%% m X n random uniform
Prob = uniform_Prob(m,n);  
%% m X n random SVD
Prob =svd_random_Prob(m,n);
%% Sparse nonsingular
density = 1/log(n*m);
rc =1/sqrt(m*n); 
Prob = sparse_matlab(m,n,density,rc);
%% Parameters
options =[];
options.ep=5*10^(-8); %norm(Prob.sol)^2/10;
options.max_time = 1;  % 10 min
options.max_iterations =  min(sum(size(Prob.A))^5, 10^10);
max_iterations = ones(4,1)*options.max_iterations ;
%% Test Distribution of method against theoretical 
pq = 0.05;
num_trials=100;
OUTPUTS ={};
method = 'Karczmarz'; %Gaussls, Gaussigp, CD, Karczmarz, LSCD 
output_stats =eval([ 'distribution_output(Prob, @iter_' method ', @boot_' method ',options,num_trials, pq );'])
output_stats.times = 1:1:length(output_stats.times);
output_theo = theoretical_output(Prob,output_stats);
OUTPUTS = [ OUTPUTS ; output_stats];
OUTPUTS = [ OUTPUTS ; output_theo];
prettyPlot_solve_system_wrapper(OUTPUTS,Prob); hold on
plot(output_stats.error_stats','--','color','red');
eval(['print -dpdf ' method '-' Prob.title '-dist.pdf' ]);
%% Testing a Single run of method against theoretical
OUTPUTS ={};
method = 'Gaussls';
[x, output_real] = eval(['solve_system(Prob,@iter_' method ', @boot_' method ',options );']);
output_real.times =  1:1:length(output_real.times);
OUTPUTS = [ OUTPUTS ; output_real];
output_theo = theoretical_output(Prob,output_real);
OUTPUTS = [ OUTPUTS ; output_theo{1}];
prettyPlot_solve_system_wrapper(OUTPUTS,Prob); 
% 
% %% All overdetermined methods against theo
% Prob.title = [Prob.title '-theo'];
% OUTPUTS ={};
% OUTPUTS= all_overdet_methods(OUTPUTS, options, Prob,4,max_iterations);
% %plot_distribution();
% % for i =1:length(OUTPUTS_mean)
% % OUTPUTS_mean{i}.times = 1:1:length(OUTPUTS_mean{i}.times);
% % end
% for i =1:length(OUTPUTS)
% OUTPUTS{i}.times = 1:1:length(OUTPUTS{i}.times);
% end
% prettyPlot_solve_system_wrapper(OUTPUTS,Prob); 
% hold on
% OUTPUTS_theos = theoretical_output(Prob,OUTPUTS);
% prettyPlot_solve_system_wrapper(OUTPUTS_theos,Prob); 
%   
% %% All P methods against theo
% Prob = uniform_Prob(m,n);  
% Prob.title = 'Gaussian-matrix';
% OUTPUTS ={};
% OUTPUTS= all_overdet_methods(OUTPUTS, options, Prob,4,max_iterations);
% for i =1:length(OUTPUTS)
% OUTPUTS{i}.times = 1:1:length(OUTPUTS{i}.times);
% end
% prettyPlot_solve_system_wrapper(OUTPUTS,Prob); 
% hold on
% %figure();
% OUTPUTS_theos = theoretical_output(Prob,OUTPUTS);
% prettyPlot_solve_system_wrapper(OUTPUTS_theos,Prob); 
