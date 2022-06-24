%---------------------------- EXPERIMENT_V -------------------------------%
%
% PLSS: (P)rojected (L)inear (S)ystems (S)olver 
%
% This experiment compares with "external" solver lsqr and two versions 
% of PLSS. Times are displayed when method converges. Otherwise
% "inf" is printed.
%
% solIdx:
%   24: PLSS    (Algorithm 1, W=I)
%   29: PLSS W  (Algorithm 1, W=1./sqrt(sum(A.*A)))
%   18: lsqr
%
% Comarison on 42 'SuiteSparse Matrix Collection' Problems.
% The problems in this experiment are small/medium overdetermined
% systems with dimensions of up to 1000 <= n  <= 10 000.
%
% Success if || Ax - b ||_2 \leq 10^{-4}
%
% This experiment is not reported in the article but run faster
% than comparisons that use random normal sketching solvers.
%
%-------------------------------------------------------------------------%
% 08/12/20, J.B., Including sequential QR method (SQR)
% 05/15/21, J.B., Including scaled solver
% 05/17/21, J.B., Modified weighted solver
% 05/07/22, J.B., Change order of solvers in experiment and preparing
%                   for release 

clc;
clear;

addpath('../ALGS');
addpath('../AUXILIARY');
addpath('../external/ssget');

% Switch warning of
warning('off','MATLAB:lsqr:tooBigTolerance');

index = ssget;

% Select rectangular matrices with column-sizes nl <= m <= nu
nl = 1000;
nu = 10000;

condit =  (index.nrows > index.ncols) & ...
          (index.ncols <= nu) & ...
          (nl <= index.ncols);

% Linear indices
ids  = find(condit);
nids = length(ids);

% Storing problem data
nprob   = 42;
names   = cell(nprob,1);
iids    = zeros(nprob,1);
nnzs    = zeros(nprob,1);
ms      = zeros(nprob,1);

% Solver data
solIdx = [24,29,18];

nsol   = length(solIdx);
exs    = zeros(nprob,nsol);
nits   = zeros(nprob,nsol);
res    = zeros(nprob,nsol);
times  = zeros(nprob,nsol);

infsS0=cell(nprob,3);

epsBase = 1e-4;
% Option(s) for PLSS
opts.print = 0; % 1
opts.tol = epsBase;

% Options(s) for 'external' solvers
M1 = [];
M2 = [];

optsS.M1 = [];
optsS.M2 = [];
optsS.itRes = [];

fprintf('********* Running EXPERIMENT IV Rect. Large ******** \n');
fprintf('*         PLSS                                       \n');
fprintf('*         Number of Solvers: %i                      \n',nsol);
fprintf('*         Number of Problems: %i                     \n',nprob);
fprintf('*         Sizes: n <= m, %i <= n <= %i               \n',nl,nu);
fprintf('*         Times in sec.                              \n');
fprintf('**************************************************** \n');
fprintf('\n');

solnames = {' PLSS    ', ' PLSS W  ',' LSQR    '};
mess = 'id \t rows \t cols \t';
for i = 1:nsol
    if i==nsol; est = '\n'; else est = '   \t';  end;
    %mess = [mess,'Sol.',num2str(i),est]; %#ok<AGROW>
    mess = [mess,solnames{i},est]; %#ok<AGROW>
end
fprintf(mess);

ts = tic;
for i = 1:nprob

    % Loading problem and dimensions
    Prob= ssget(ids(i));
    A   = Prob.A;

    [m,n] = size(A);    
    ids(i)   = Prob.id;
    names{i} = Prob.name;
    nnzs(i) = nnz(A);
    ms(i) = m;

    % Solution, tolerances and initialization
    x0_ = ones(n,1);
    x0_(1) = 10;
    b = A*x0_;
    nb = norm(b);
    x0 = zeros(n,1);

    resEps = epsBase/nb;

    MAXIT = n;
    opts.maxiter = MAXIT;
    opts.tol     = epsBase;

    optsS.tolRel  = resEps;
    optsS.maxit = MAXIT;    
    optsS.optsP = opts;
    
    % Store true solution for possible comparisons
    optsS.x=x0_;

    % Run solvers and store outcomes
    for j = 1:nsol

        outs = runLSolver(A,b,x0,optsS,solIdx(j));
        if outs.res < epsBase; exs(i,j) = 1; else exs(i,j) = 0; end
        nits(i,j)   = outs.niter;
        res(i,j)    = outs.res;
        times(i,j)  = outs.time;
        
    end

    messb = '%i\t%i\t %i ';
    messm = repmat('\t %.2e     ',1,(nsol));
    mess = [messb,messm,'\n'];
    fprintf(mess,i,m,n,times(i,1:nsol)./exs(i,1:nsol));

end
te = toc(ts);

% Switch warning message on
warning('on','MATLAB:lsqr:tooBigTolerance');
