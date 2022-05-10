%%Test methods for overdetermined systems including
% Randomized Kaczmarz
% Coordinate Descent for Least Squares
% Gaussian Kaczmarz
% Gaussian for Least Squares
% Robert Gower and Peter Richtarik -- 18 June 2015
%% Select Dimensions
clear all
m = 200;  % number of rows
n =100;    %select dimension
%% Badly conditioned positive definites  % uncomment to test a dense matrix
% Prob = uniform_Prob(m,n);
%% Sparse nonsingular   % uncomment to test a sparse matrix
% density = 1/log(n*m); rc = 1/sqrt(m*n); 
% Prob = sparse_matlab(m,n,density,rc);
%% MatrixMarket Format
clear Prob
Prob.title= 'illc1033';
[Prob.A,rows,cols,entries,rep,field,symm] = mmread([Prob.title '.mtx']); 
[Prob.b,~,~,~,~,~,~] = mmread([Prob.title '_rhs1.mtx']); 
Prob.sol = Prob.A\Prob.b;
%% Select method Parameters
if(~isfield(Prob,'sol'))
    Prob.sol = Prob.A\Prob.b;
end
options =[];
options.ep=1*10^(-4);
options.max_time = 5*60;  % 300s
options.max_iterations =  10^6;  % maximum of 10^6 iterations
options.metric = Prob.A'*Prob.A;  % this forces all methods to use the same error measurement, comment this, the methods will select their own pd matrices for measuring error
OUTPUTS ={};
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
%figure('visible','on')
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)
%% Plot all methods and 95% and 5% as shaded regions
 %plot_distribution();



