%---------------------------- EXPERIMENT_II ------------------------------%
%
% PLSS: (P)rojected (L)inear (S)ystems (S)olver 
%
% This experiment compares with "external" solver lsqr, a random normal
% projection method and two versions of PLSS. 
% Times are displayed when method converges. Otherwise "inf" is printed.
%
% solIdx:
%   30: Random normal sketch
%   24: PLSS    (Algorithm 1, W=I)
%   28: PLSS W  (Algorithm 1, W=1./sqrt(sum(A.*A)))
%   18: lsqr
%
% Comarison on 51 'SuiteSparse Matrix Collection' Problems.
% The problems in this experiment are large overdetermined
% systems with dimensions of up to 10 000 <= n <= Inf.
%
% Success if || Ax - b ||_2 \leq 10^{-2}
%
% This experiment is Table 2 in the corresponding article
%
% Note: This experiment is relatively time consuming
%
%-------------------------------------------------------------------------%
% 05/17/21, J.B., Including weighted modified solver
% 12/07/21, J.B., Including W=I solver
% 12/08/21, J.B., Including sparse random solver, omitting very large
%                   instances for random solver
% 05/07/22, J.B., Preparing files for release

clc;
clear;

exName ='II';

addpath('../ALGS');
addpath('../AUXILIARY');
addpath('../../../../external/ssget');

% Switch warning of
warning('off','MATLAB:lsqr:tooBigTolerance');

index = ssget;

% Select rectangular matrices with column-sizes nl <= n <= nu
nl = 10000;
nu = Inf;

condit =  (index.nrows > index.ncols) & ...
          (index.ncols <= nu) & ...
          (nl <= index.ncols);

% Linear indices
ids  = find(condit);
nids = length(ids);

% Storing problem data
nprob   = 51; % Last problem very large and uses much time
names   = cell(nprob,1);
iids    = zeros(nprob,1);
nnzs    = zeros(nprob,1);
ms      = zeros(nprob,1);
ns      = zeros(nprob,1);
nxs      = zeros(nprob,1);

% Solver data
solIdx = [30,24,28,18];% 27,
nsol   = length(solIdx);
exs    = zeros(nprob,nsol);
nits   = zeros(nprob,nsol);
res    = zeros(nprob,nsol);
times  = zeros(nprob,nsol);

infsS0=cell(nprob,3);

epsBase = 1e-2;
% Option(s) for PLSS
opts.print = 0;
opts.tol = epsBase;

% Flags for additional storage
opts.storeSolRes=1;
opts.storeRelRes=1;

%Additioanal storage
relRes  = zeros(nprob,nsol);
solRes  = zeros(nprob,nsol);

% Options(s) for 'external' solvers
M1 = [];
M2 = [];

optsS.M1 = [];
optsS.M2 = [];
optsS.itRes = [];

fprintf('********* Running EXPERIMENT II Rect. Large ******** \n');
fprintf('*         PLSS                                       \n');
fprintf('*         Number of Solvers: %i                      \n',nsol);
fprintf('*         Number of Problems: %i                     \n',nprob);
fprintf('*         Sizes: n <= m, %i <= n <= %i               \n',nl,nu);
fprintf('*         Times in sec.                              \n');
fprintf('**************************************************** \n');
fprintf('\n');

solnames = {' Rnd. Prj',' PLSS    ', ' PLSS W  ',' LSQR    '};
mess = 'id \t rows \t cols \t';
for i = 1:nsol
    if i==nsol; est = '\n'; else est = '   \t';  end;
    %mess = [mess,'Sol.',num2str(i),est]; %#ok<AGROW>
    mess = [mess,solnames{i},est]; %#ok<AGROW>
end
fprintf(mess);

ts = tic;
%41 is very big, 43 and 52 are huge
for i = 1:nprob % i=20:20%

    Prob= ssget(ids(i));
    A   = Prob.A;

    [m,n] = size(A);
    ids(i)   = Prob.id;
    names{i} = Prob.name;
    nnzs(i) = nnz(A);
    ms(i) = m;
    ns(i) = n;

    % Problem
    x0_ = ones(n,1);
    x0_(1) = 10;
    b = A*x0_;
    nb = norm(b);
    x0 = zeros(n,1);

    nxs(i) = norm(x0_);
    
    resEps = epsBase/max(nb,1);

    optsS.x = x0_;
    
    % Cap the maximum iterations for large problems
    MAXIT = min(n,500); % 500
    
    opts.maxiter = MAXIT; %500;
    opts.tol     = epsBase;

    optsS.tolRel  = resEps;
    optsS.maxit = MAXIT; %100;
    optsS.optsP = opts;
    
    optsS.storeSolRes=opts.storeSolRes;
    optsS.storeRelRes=opts.storeRelRes;
    
    % For very large problems use smaller projection
    if n>1e5
        optsS.r=5;
    else
        optsS.r=min(50,0.1*n);
    end
    
    for j = 1:nsol

        outs = runLSolver(A,b,x0,optsS,solIdx(j));
        if outs.res < epsBase; exs(i,j) = 1; else exs(i,j) = 0; end
        nits(i,j)   = outs.niter;
        res(i,j)    = outs.res;
        times(i,j)  = outs.time;
       
        if isfield(opts,'storeRelRes')
            if opts.storeRelRes==1;
                relRes(i,j)=outs.relRes;
            end
            
        end
                
        if isfield(opts,'storeSolRes')
            if opts.storeSolRes==1
                solRes(i,j)=outs.solRes;
            end
        end
        
    end

    messb = '%i\t%i\t %i ';
    messm = repmat('\t %.2e     ',1,(nsol));
    mess = [messb,messm,'\n'];
    fprintf(mess,i,m,n,times(i,1:nsol)./exs(i,1:nsol));

end
te = toc(ts);

savepath = fullfile('./data');
save([savepath,'/EXPERIMENT_',exName],'names','iids','nnzs','ms','ns',...
    'exs','nits','res','times','solRes','relRes'); % ,'solRes','relRes'

