%---------------------------- EXPERIMENT_III -----------------------------%
%
% PLSS: (P)rojected (L)inear (S)ystems (S)olver 
%
% This experiment compares with "external" solver craig, a random normal
% projection method and two versions of PLSS. 
% Times are displayed when method converges. Otherwise "inf" is printed.
%
% solIdx:
%   30: Radom normal sketch
%   24: PLSS    (Algorithm 1, W=I)
%   28: PLSS W  (Algorithm 1, W=1./sqrt(1+sum(A.*A)))
%   34: craig
%
% Comarison on 42 'SuiteSparse Matrix Collection' Problems.
% The problems in this experiment are small/medium overdetermined
% systems with dimensions of up to 1000 <= n  <= 10 000.
%
% Success if || Ax - b ||_2 \leq 10^{-4}
%
% This experiment is Table 3 in the corresponding article
%
%-------------------------------------------------------------------------%
% 08/12/20, J.B., Including sequential QR method (SQR)
% 05/15/21, J.B., Including scaled solver
% 05/17/21, J.B., Modified weighted solver
% 05/19/21, J.B., Including randomized solver in comparisons
% 12/03/21, J.B., Including PLSS W=I in comparison
% 12/06/21, J.B., Storing relative residuals and distance to solution
% 05/07/22, J.B., Preparing for release
% 06/24/22, J.B., Further preparations

clc;
clear;

exName ='III';

addpath('../ALGS');
addpath('../AUXILIARY');
addpath('../external/ssget');
addpath('../external/craigSOL');

index = ssget;

% Select rectangular matrices with column-sizes nl <= n <= nu

nl = 1000;
nu = 10000;

condit =  (index.nrows < index.ncols) & ...
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
ns      = zeros(nprob,1);

% Solver data
solIdx = [30,24,28,34];
nsol   = length(solIdx);
exs    = zeros(nprob,nsol);
nits   = zeros(nprob,nsol);
res    = zeros(nprob,nsol);
times  = zeros(nprob,nsol);

infsS0=cell(nprob,3);

epsBase = 1e-4; % 1e-2, 1e-4
% Option(s) for PLSS
opts.print = 0; % 1
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

fprintf('********* Running EXPERIMENT III Underdet. ******** \n');
fprintf('*         PLSS                                      \n');
fprintf('*         Number of Solvers: %i                     \n',nsol);
fprintf('*         Number of Problems: %i                    \n',nprob);
fprintf('*         Sizes: n <= m, %i <= n <= %i              \n',nl,nu);
fprintf('*         Times in sec.                             \n');
fprintf('*************************************************** \n');
fprintf('\n');

solnames = {' Rnd. Prj',' PLSS    ', ' PLSS W  ',' CRAIG    '};
mess = 'id \t rows \t cols \t';
for i = 1:nsol
    if i==nsol; est = '\n'; else est = '   \t';  end;
    % mess = [mess,'Sol.',num2str(i),est]; %#ok<AGROW>
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
    ns(i) = n;

    % Solution, tolerances and initialization
    x0_ = ones(n,1);
    x0_(1) = 10;
    b = A*x0_;
    
    nb = max(norm(b),1);    
    x0 = zeros(n,1);
    resEps = epsBase/nb;
        
    MAXIT = n + 1500;       
    opts.maxiter = MAXIT;
    opts.tol = epsBase; % epsBase
    
    % Set r (number of columns in the random solver)
    r = 10; % 10, 50
    opts.r = r;
    
    optsS.tolRel  = resEps;
    optsS.maxit = MAXIT;
    optsS.optsP = opts;
    nproj = 1;
    nintsol = 1;
    
    % Parameters for Craig's method
    optsS.atol   = 0; %1e-10; %opts.tol;
    optsS.btol   = resEps; %opts.tol;
    optsS.conlim = 1.0e+14;    
    optsS.show   = 0;   
    
    % True solution (for later comparison)
    optsS.x=x0_;
%     %opts.storeResSoln=1;
%     %opts.storeRelRes=1;
% 
%     % Adjust limited memory parameters
%     memRat = 1/100;
%     mem = 25;%min(floor(memRat*n),100);
%     nRst = 300;%2*mem;
% 
%     optsS.optsP.rst = nRst;
%     optsS.optsP.mem = mem;
%     optsS.itRes = nRst;
%     
    optsS.storeSolRes=opts.storeSolRes;
    optsS.storeRelRes=opts.storeRelRes;

    for j = 1:nsol

%         if solIdx(j) == 21
%             %optsS.optsP.maxiter = ;
%             optsS.optsP.rst = MAXIT;
%             optsS.optsP.mem = MAXIT;
%         end
        
        outs = runLSolver(A,b,x0,optsS,solIdx(j));
        if outs.res < epsBase; exs(i,j) = 1; else exs(i,j) = 0; end % epsBase
        nits(i,j)   = outs.niter;
        res(i,j)    = outs.res;
        times(i,j)  = outs.time;
%         if solIdx(j) == 12 || solIdx(j) == 13
%             infsS0{i,nintsol} = outs.info1;
%             nintsol = nintsol + 1;
%         end
        
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
    'exs','nits','res','times','solRes','relRes');

