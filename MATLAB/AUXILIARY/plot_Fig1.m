%------------------------------- plot_Fig1.m -----------------------------%
%
% Script to replot experiment data (called Fig.1)
%
%-------------------------------------------------------------------------%
% 05/19/21, J.B., Preparation of file for using data
% 05/08/22, J.B., Preparation for release
clc;
clear;

figpath = fullfile(pwd,'..','EXPERIMENTS/figs/');
datapath = fullfile(pwd,'..','EXPERIMENTS/data/');

indAlg = [1 2 3 4];
np = 42;

nsol = 4;

% Plotting labels, and additional parameters
leg={   'Rand. Proj.',...
        'PLSS',...
        'PLSS W',...
        'LSQR'};

%mleg = cell(nms);
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

% Initial line types
types.colors    = ['r' 'b' 'b' 'k']; %'k' 'y'
types.lines     = {'-', '-.', '-','-','-','-.','-'}; %'-',   '-'
leglocation = 'SouthEast';
isLog = 1;

exName ='I';
name_EXPERIMENT = ['EXPERIMENT_',exName];
data = load([datapath,name_EXPERIMENT],'exs','nits','times');
exs = data.exs;
times = data.times;

% Specification of tick-marks for extended perf. profiles
ticks   = -8; 
ticke   = 24; 
XTick   = 2.^(ticks:ticke);
XLim    = [XTick(1),XTick(end)];
leglocation1 = 'SouthEast';
legFontSize1 = 10;
markerSize = 0;
lineWidths = [1.8, 1.8, 1.8 1.8];

% Extended perf. profile
perf_ext_fnc(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types,...
    leglocation1,XTick,XLim,legFontSize1,markerSize,lineWidths);

% Axis annotation
ax = gca;
ax.XTick = [2^(-8), 2^(-4), 1, 2^4, 2^8, 2^(12), 2^16, 2^20, 2^(24)];
XTickLabel = {'2^{-8}', '2^{-4}', '1', '2^{4}', '2^{8}', '2^{12}', '2^{16}' , '2^{20}' '2^{24}'};
set(ax,'XTickLabel',XTickLabel);

% On first plot "grid on" does not produce desired result,
% thus toggle once more
grid on; grid off;
grid on;
title('$\textnormal{Times}$','Interpreter','latex');

% Modified printing
fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];
figname = ['Fig1','.pdf'];

print(fig,'-dpdf',fullfile(figpath,figname));

close ALL;
