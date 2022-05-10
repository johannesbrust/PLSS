%%Test methods for overdetermined systems including
% Randomized Kaczmarz
% Coordinate Descent for Least Squares
% Gaussian Kaczmarz
% Gaussian for Least Squares
% Robert Gower and Peter Richtarik -- 18 June 2015

%------------------------- Modifications ---------------------------------%
%
% 07/01/20, J.B., Including PLSSr (full memory)
% 09/06/20, T.M., Include PLSSRandom (fixed sample size AND increasing sample size)
%                 Include modif of the 2 codes to use 'calculate_error'
% 11/12/20, T.M., add cleaning of the max time, as the boot step can take several seconds.
% 12/07/20, T.M., uniform stop
%
%-------------------------------------------------------------------------%
%% Select Dimensions
clear all
m = 500000;  % number of rows
n =200;    %select dimension
%% Badly conditioned positive definites  % uncomment to test a dense matrix
 Prob = uniform_Prob(m,n);nn=n;
%% Sparse nonsingular   % uncomment to test a sparse matrix
% density = 1/log(n*m); rc = 1/sqrt(m*n); 
% Prob = sparse_matlab(m,n,density,rc);
%% MatrixMarket Format
%{
clear Prob
Prob.title= 'illc1033';
[Prob.A,rows,cols,entries,rep,field,symm] = mmread([Prob.title '.mtx']); 
[m,nn] = size(Prob.A);
%
[Prob.b,~,~,~,~,~,~] = mmread([Prob.title '_rhs1.mtx']); 
Prob.sol = Prob.A\Prob.b;
%}

%%%%%%%%%%%% Tangi %%%%%%%%%%%%%%%
% Create a consistent system
    x0_ = ones(nn,1);
    x0_(1) = 10;
    Prob.b = Prob.A*x0_;
    Prob.sol = x0_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% Select method Parameters
if(~isfield(Prob,'sol'))
    Prob.sol = Prob.A\Prob.b;
end
options =[];
options.ep=1*10^(-6);
options.max_time = 1;%5*60;  % 300s
options.max_iterations =  10^6;nn;%10^6;  % maximum of 10^6 iterations
options.metric = eye(size(Prob.A,2)); %Prob.A'*Prob.A;  % this forces all methods to use the same error measurement, comment this, the methods will select their own pd matrices for measuring error
OUTPUTS ={};
%% Gaussian Least squares
[x, output_gaussls] = solve_system(Prob,@iter_Gaussls, @boot_Gaussls,options );
output_gaussls.times = min(output_gaussls.times, options.max_time);
%OUTPUTS = [ OUTPUTS ; output_gaussls];
%% LS Coordinate Descent
[x, output_lscd] = solve_system(Prob,@iter_LSCD, @boot_LSCD,options );
output_lscd.times = min(output_lscd.times, options.max_time);
%OUTPUTS = [ OUTPUTS ; output_lscd];
%% Gaussian Iterative Projection IGP
[x, output_gaussigp] = solve_system(Prob,@iter_Gaussigp, @boot_Gaussigp,options );
output_gaussigp.times = min(output_gaussigp.times, options.max_time);
%OUTPUTS = [ OUTPUTS ; output_gaussigp];
%% Randomized Karczmarz
[x, output_kz] = solve_system(Prob,@iter_Karczmarz2, @boot_Karczmarz2,options );
output_kz.times = min(output_kz.times, options.max_time);
OUTPUTS = [ OUTPUTS ; output_kz];

%% PLSS
opts.tol = 1e-15;
opts.maxtime = options.max_time;%5*60;  % 300s
opts.maxiter = 10^6;nn;%10^6;  % maximum of 10^6 iterations
opts.metric = eye(size(Prob.A,2));
opts.projRes = false;
opts.store      = 1;
opts.errorScale = 100;
opts.ep         = 1*10^(-6);
%opts.print =    1;
x0              = zeros(nn,1);
opts.zeroStart  = 1; %zero is the starting point
%{
[xkp,nckp,outs] = plss_1_scaledRLS(x0,Prob.A,Prob.b,opts);
output_plss.errors  = outs.errs';
output_plss.times   = outs.times';
output_plss.name    = 'PLSS';
OUTPUTS = [ OUTPUTS ; output_plss];
%}
%opts.rst  = 250;
%opts.mem  = 10;
opts.sampleSize = 1; %10% of nn
[ xkr, nckr, outsr ] = plss_1_rand(x0,Prob.A,Prob.b, opts, Prob );
output_plssr.errors  = outsr.errs';
output_plssr.times   = outsr.times';
output_plssr.name    = 'PLSSR-1';
OUTPUTS = [ OUTPUTS ; output_plssr];
%{
opts.sampleSize = 10; %10% of nn
[ xkr, nckr, outsr ] = plss_1_rand(x0,Prob.A,Prob.b, opts, Prob );
output_plssr.errors  = outsr.errs';
output_plssr.times   = outsr.times';
output_plssr.name    = 'PLSSR-10';
opts.sampleSize = 100; %10% of nn
[ xkr, nckr, outsr ] = plss_1_rand(x0,Prob.A,Prob.b, opts, Prob );
output_plssr.errors  = outsr.errs';
output_plssr.times   = outsr.times';
output_plssr.name    = 'PLSSR-100';
%OUTPUTS = [ OUTPUTS ; output_plssr];
[ xkri, nckri, outsri ] = plss_1_rand_increase(x0,Prob.A,Prob.b, opts, Prob );
output_plssri.errors  = outsri.errs';
output_plssri.times   = outsri.times';
output_plssri.name    = 'PLSSRi';
OUTPUTS = [ OUTPUTS ; output_plssri];
%}

%% plotting
%figure('visible','on')
prettyPlot_solve_system_wrapper(OUTPUTS,Prob)
%% Plot all methods and 95% and 5% as shaded regions
 %plot_distribution();



