%------------------------ example_2 --------------------------------------%
% 
% Example of applying PLSS (Algorithm  1, W=I)
% to an ill-conditioned nonsymmetric matrix
% The matrix is build-in in Matlab
%
% n = 100;
% A = sampling(n);
%
% Random normal projections are included for comparison (when size of
% sketch varies)
%
% The solution is x = ones(n,1);
%-------------------------------------------------------------------------%
% 05/08/22, J.B., Preparation for release

clc;
clear;
close all;

addpath('../ALGS');
figpath = fullfile(pwd,'..','EXPERIMENTS/figs/');

%---------- Initialize figure size
width =4; height =3;

figure('Units','inches',...
        'Position',[0 0 width height],...
        'PaperPositionMode','auto');
%----------    

% Setup problem
n       = 100;
A       = gallery('sampling',n);

[m,~]   = size(A);

b       = full(sum(A,2));
nb      = norm(b);
x0      = zeros(n,1);

nrand   = 3;
nsol    = nrand+1;
maxit   = n;
tol     = 1e-5;
rs      = [10,20,40]; 

fprintf('********* Running Example 2 *********************** \n');
fprintf('*         PLSS Solvers                              \n');
fprintf('*         Number of Solvers: %i                     \n',nsol);
fprintf('*         Sizes: m=%i, n=%i                         \n',m,n);
fprintf('*         Tol.= %1.6f                               \n',tol);
fprintf('*************************************************** \n');
fprintf('\n');

% Solver options
opts.tol = tol;
opts.print = 0;
opts.store = 1;
opts.maxiter = maxit;

% Call PLSS solver
[xk,nck,outs] = plss_R(x0,A,b,opts);

% Plot outcomes
msp = 12;
nres = length(outs.errs);
hold on;
semilogy(0:nres-1,outs.errs,'-','LineWidth',2,'color','cyan');

ftnsize     = 15;
ftnsizeL    = 12;
lw = 1;
names = cell(nsol,1);
name1 = '$\textnormal{PLSS}$';
names{1} = name1;
for i=1:nrand
    
    % Call Random normal solver
    r = rs(i);
    opts.r = r;
    [xkRND,nckRND,outsRND] = plss_RANDN(x0,A,b,opts);    
    nres1 = length(outsRND.errs);
    
    % Plot next outcome
    semilogy(0:nres1-1,outsRND.errs,'-','LineWidth',lw);
    name2       = sprintf('$\\textnormal{Rand. Proj. } (r=%i)$',r);
    names{i+1}  = name2;
end
    
legend(names,...
    'FontSize',ftnsizeL,'Location','North','Interpreter','latex');
xlabel('$k$','FontSize',ftnsize,'Interpreter','latex');
ylabel('$\|\mathbf{r}_k\|_2$',...
    'FontSize',ftnsize,'Interpreter','latex');
box on;

% Save file
% On first plot "grid on" does not produce desired result,
% thus toggle once more
%grid on; grid off;
%grid on;
title('$\textnormal{Projection Solvers}$','Interpreter','latex',...
    'FontSize',ftnsize);

set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 0.5         );

% ax = gca;
% ax.YTick = [0 500 1000];
% YTickLabel = {'0', '500', '1000'};
% set(ax,'YTickLabel',YTickLabel);
% Modified printing
fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];
figname = ['Fig2','.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));
