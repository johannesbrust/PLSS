function prettyPlot_solve_system_wrapper_figpath(OUTPUTS,Prob,figpath)
lO = length(OUTPUTS);
ERRORS = {};
TIMES = {};
legendStr = {};
for i =1:lO
    lent(i) = length(OUTPUTS{i}.times);
    if(lent(i)>1000)
        grid =1:floor(lent(i)/1000):lent(i);
        lent(i)=length(grid);
    else
        grid = 1:1:lent(i);
    end
    ERRORS = [ERRORS ; OUTPUTS{i}.errors(grid)];
    TIMES = [TIMES; OUTPUTS{i}.times(grid)];
    legendStr = [legendStr ; OUTPUTS{i}.name];
end
legendStr
markers={'o','+','^','*','s','d','v','+','<','>','x','h','.',...
'^','*','o','p','+','<','h','.','>','x','s','d','v',...
'o','p','+','*','s','d','v','+','<','>','x','h','.'};
options.title = Prob.title;
options.markersize =8;
options.linewidth =1.5;
options.markers= markers(1:lO);
options.colors = lines(lO);
options.xlabel = '$\textnormal{time (s)}$';
options.ylabel ='$\| \mathbf{r}_k \|_2$';
options.ylim = [-0.01 10];
options.xlim = [-0.01 1.01];
options.logScale = 2;
markerspace = [ceil(lent/7)' ceil(rand(lO,1).*(lent./7)')];
options.markerSpacing = markerspace;
lineStyles= {':','--','-','--','-','-','-','-'};
options.lineStyles = lineStyles(1:lO);
options.legend = legendStr;
options.legendLoc = 'SouthEast';
%h= figure();
% x0=10;y0=10;width=500;height=350;
% set(gca, 'Position', [x0 y0 width height])
options.title(ismember(options.title,'_ ,.:;!/\')) = '-';
prettyPlot(TIMES,ERRORS, options);

% Adjust axis after plot
%axis tight

options.ylim = [-0.01 10];
options.xlim = [-0.01 1.01];

options.title = [ 'figures/' options.title ]
%set(gca,'defaultTextFontName', 'Arial')
set(gca,'fontsize',20)
xlhand = get(gca,'xlabel');  set(xlhand,'fontsize',20) ;
xlhand = get(gca,'ylabel');  set(xlhand,'fontsize',20) ;
xlhand = get(gca,'title');  set(xlhand,'fontsize',15) ;
ch = get(gcf,'children'); set(ch(1), 'fontsize',13);
if (~isempty(findstr('interface',pwd)))
    cd ..
end

%print -dpdf
print(figpath,'-dpdf');

% if (isfield(options,'logScale'))
%     eval(['print -dpdf ' options.title '.pdf' ]);
% else
%     eval(['print -dpdf ' options.title '-not-log' '.pdf' ]);
% end
%print -dpdf prettyPlot.jpg
end