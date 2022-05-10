%----------------------- EXPERIMENT_III ----------------------------------%
%
% Comparison of PLSS(r) to methods from Gower, Richtarik (GR), (2015)
%
% This test uses LIBSVM problems in order to do "Ridge Regressions"
%
% Extended experiment with randomized methods
%
% Solvers in the comparison are
%
%   Gauss pd (GR)
%   CD pd (GR)
%   Block CD-pd (GR)
%   Block Gauss (GR)
%   PLSS W (Algorithm 1, W=1./sqrt(sum(A.*A)))
%   PLSS (Algorithm 1, W=I)
%   PLSS KZ (Kaczmarz like method)
%
% Comarison on 8 LIBSVM Problems.
% The problems and setup in in this experiment are the same from GR.
%
% Success if || Ax - b ||_2 \leq 10^{-4}
%
% This experiment corresponds to Figure 1 in the article of
% J.J. Brust and M.A Saunders (2022)
%-------------------------------------------------------------------------%
% 07/01/20, J.B.
% 07/02/20, J.B., modifications for comparisons (e.g., error measure in gr
% methods.
% 07/07/20, J.B., Running GR solvers for longer times
% 07/20/20, T.M., add random PLSS solver
% 12/08/21, J.B., Test run of experiment and using plss_RW1
% 01/11/22, J.B., Including randomized methods
% 04/28/22, J.B., Including extended Kaczmark method
% 05/03/22, J.B., Randomized PLSS Kaczmark method
% 05/07/22, J.B., Preparing for release

clc;
clear;
close all;

addpath(genpath('../ALGS'));

% Specify problem path
myp = strsplit(pwd,'/');
myp{end} = 'external';
extPath = fullfile(myp{:});

libSVMpath = ['/',extPath,'/','LIBSVM']; % /libsvm/matlab
RandLinLab = ['/',extPath,'/','RandomLinearLab'];
addpath(genpath(libSVMpath));
addpath(genpath(RandLinLab));

figpath ='figs/';

probPath = [libSVMpath,'/data/'];

problist = {
    'a6a';...
    'a7a';...
    'a8a';...
    'a9a';...
    'aloi.scale';...
    'protein';...
    'covtype.libsvm.binary.scale';...
    'SUSY'};

nprobs = length(problist);

timesr = zeros(nprobs,1);
timesp = zeros(nprobs,1);

% Randomized solvers options
options =[];
options.ep=1*10^(-4);
options.max_time = 1;%1%5*60;  % 300s

nsol = 8;

fprintf('********* Running EXPERIMENT III (Randomized) ******** \n');
fprintf('*         PLSS                                         \n');
fprintf('*         Number of Solvers: %i                        \n',nsol);
fprintf('*         Number of Problems: %i                       \n',nprobs);
fprintf('*         Conv.: norm(rk) <= %1.1e                     \n',options.ep);
fprintf('*         Max time = %1.f sec.                         \n',options.max_time);
fprintf('****************************************************** \n');
fprintf('\n');


%% Read LIBSVM problems.
for i = 1:nprobs % 3:3%
    
    file    = problist{i};
    
    % Reading problem
    tr = tic;
    if strcmp(file,'SUSY')==1
        SUSY_AB = load('./data/SUSY_AB');
        AA = SUSY_AB.AA;
        Ab = SUSY_AB.Ab;
        m = SUSY_AB.m;
        n = SUSY_AB.n;
    else
        [b,A]   = libsvmread([probPath,file]);
        timesr(i) = toc(tr);
        
        % Forming problem
        [m,n] = size(A);
        In = eye(n);
        tp = tic;
        AA = (A'*A)+In;
        Ab = A'*b;
        timesp(i) = toc(tp);
        
    end
        
    clear Prob;
    %close all;
    Prob.A = AA;
    Prob.b = Ab;
    Prob.sol = Prob.A\Prob.b;
    Prob.title = file;
    nb = norm(Ab);
    
    %%----------
    % Solvers from RandomLinearLab
    %
    % Gausspd
    % CD
    % BCD
    % GaussBpd
    %
    %-----------
    options.errorScale = nb^2;
    options.errorSqrt = 1;
    % Calling solvers
    options.max_iterations =  10^6;%n;%10^6;  % maximum of 10^6 iterations
    
    % This forces all methods to use the same error measurement, 
    % comment this, the methods will select their own pd matrices for measuring error
    options.metric = (Prob.A'*Prob.A);% In; %Prob.A'*Prob.A;  
    OUTPUTS ={};
    [~, output_gausspd] = solve_system(Prob,@iter_Gausspd, @boot_Gausspd,options );
    OUTPUTS = [ OUTPUTS ; output_gausspd]; %#ok<*AGROW>
    %% Randomized Coordinate Descent for positive definite
    [~, output_cd] = solve_system(Prob,@iter_CD, @boot_CD,options );
    OUTPUTS = [ OUTPUTS ; output_cd];
    %% Block Coordinate Descent for positive definite
    [~, output_Bcd] = solve_system(Prob,@iter_BCD, @boot_BCD,options );
    OUTPUTS = [ OUTPUTS ; output_Bcd];
    %% %% BLOCK METHODS
    %% Block Gaussian positive definite
    [~, output_gausspd] = solve_system(Prob,@iter_GaussBpd, @boot_GaussBpd,options );
    output_gausspd.name = '$\textnormal{Block Gauss}$';
    OUTPUTS = [ OUTPUTS ; output_gausspd];
    
    %%----------
    % Solvers from PLSS
    %
    % PLSS W 
    % PLSS
    % PLSS KZ
    %
    %-----------
    %% PLSS W
    opts.projRes = false;
    opts.maxiter    = max(m,n);
    opts.store      = 1;
    opts.tol        = options.ep;
    x0              = zeros(n,1);
    opts.print = 0;
    [~,~,outs] = plss_RW1(x0,Prob.A,Prob.b,opts);
    le = length(outs.errs);
    stp = 1;
    if le > 1e3 && le > n
        stp = 50;
    elseif le > n
        stp = 10;
    end
    output_plss.errors  = outs.errs';
    output_plss.errors  = output_plss.errors([1:stp:le,le]);
    output_plss.times   = outs.times';
    output_plss.times   = output_plss.times([1:stp:le,le]);
    output_plss.name    = 'PLSS W';
    OUTPUTS = [ OUTPUTS ; output_plss];
    
    %% PLSS
    opts.projRes = false;
    opts.maxiter    = max(m,n);
    opts.store      = 1;
    opts.tol        = options.ep;
    x0              = zeros(n,1);
    opts.print = 0;
    [~,~,outs] = plss_R(x0,Prob.A,Prob.b,opts);
    
    le = length(outs.errs);
    stp = 1;
    if le > 1e3 && le > n
        stp = 50;
    elseif le > n
        stp = 10;
    end
    output_plss.errors  = outs.errs';
    output_plss.errors  = output_plss.errors([1:stp:le,le]);
    output_plss.times   = outs.times';
    output_plss.times   = output_plss.times([1:stp:le,le]);
    output_plss.name    = 'PLSS';
    OUTPUTS = [ OUTPUTS ; output_plss];
           
    %% PLSS KZ
    opts.projRes = false;
    opts.maxiter    = min(m,n+50);
    opts.store      = 1;
    opts.tol        = options.ep;
    x0              = zeros(n,1);
    opts.print = 0;
    
    [~,~,outs] = plss_KZ_EX(x0,Prob.A,Prob.b,opts);
    
    le = length(outs.errs);
    stp = 1;
    if le > 1e3 && le > n
        stp = 50;
    elseif le > n
        stp = 10;
    end
    output_plss.errors  = outs.errs';
    output_plss.errors  = output_plss.errors([1:stp:le,le]);
    output_plss.times   = outs.times';
    output_plss.times   = output_plss.times([1:stp:le,le]);    
    output_plss.name    = 'PLSS KZ';
    OUTPUTS = [ OUTPUTS ; output_plss];
    
    figpn = [figpath,file,'.pdf'];
    figure;
    prettyPlot_solve_system_wrapper_figpath(OUTPUTS,Prob,figpn);
    
end


