function perf_ext_fnc(ex,T,legstr, varargin)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004
%
% Modified by Spartak Zikrin, November 2015

% 12/11/16, J.B. extension for plotting on beamer slides.
% 07/12/17, J.B. extension for passing line types.
% 05/02/19, J.B. extension for passing legend location.
% 09/03/19, J.B. extension for passing tau end values.
% 01/22/20, J.B., Extended performance profile based on Leyffer et al., '12

[np,ns] = size(T);


if nargin < 3
    leg=1:ns;
    leg=leg';
    legstr=num2str(leg);
end

isbeam = 0;
if nargin > 3
    
    isbeam = varargin{1};
    
end

colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   
lines   = {'-' '--' '-.'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];

leglocation = 'SouthEast';

hastypes = 0;
if nargin > 4
   
    hastypes = 1;
    types   = varargin{2};
    colors  = types.colors;
    lines   = types.lines;
    markers = types.markers;
    
end
if nargin > 5
    
    leglocation = varargin{3};
    
end

hasLimits = false;
if nargin > 6
   XTick    = varargin{4};
   XLim     = varargin{5};
   hasLimits = true;
end

if ns > 3
    legFontSize = 8;
else
    legFontSize = 14;
end
if nargin > 7
    
    legFontSize = varargin{6};
    
end

markerSize = 6;
useMarkers = 0;
if nargin > 8
    
    markerSize = varargin{7};
    if markerSize > 0
        useMarkers = 1;
    end
    
end

hasLineWidth = 0;
if nargin > 9
    
    lineWidths = varargin{8};
    hasLineWidth = 1;
    
end


T(ex~=1)=NaN;

% Minimal performance per solver
%minperf = min(T,[],2);
sidx    = 1:ns;

% Compute ratios and divide by smallest element in each row.
% Modification for extended performance profile. Take minimum over
% all other solver results to compute ratio.
r = zeros(np,ns);
for p = 1: np
    
    for j = 1:ns
        
        tmpmin = min(T(p,sidx(sidx~=j)));
        
        if ~isnan(tmpmin) 
        
            r(p,j) = T(p,j)/tmpmin;
            
        else
           
            r(p,j) = T(p,j)/T(p,j);
            
%             if ~isnan(T(p,j))
%                 
%                 r(p,j) = 1.0;
%                 
%             end            
            
        end
        
        %r(p,:) = T(p,:)/minperf(p);
    end
end

if hasLimits == true
    max_ratio = XLim(2);
else
    max_ratio = max(max(r));
end

%max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.
r(isnan(r)) = 2*max_ratio;
r = sort(r);

% Plot stair graphs with markers.
% Set figure properties
if isbeam == 1
    x0 = 0; y0 = 0;
    width =4; height =3;

    figure('Units','inches',...
        'Position',[x0 y0 width height],...
        'PaperPositionMode','auto');
else
    figure;
end

hold on;

% legend(legstr,'Location','SouthEast');
ymax=1;
ymin=1;
for s = 1: ns
    [xs,ys] = stairs(r(:,s),[1:np]/np);
    if ys(end)==1
         ys(end+1)=1;
         ymax=1.01;
         xs(end+1)=1.05*max_ratio;
    end
    
    if hastypes == 0
        sl = mod(s-1,3) + 1; 
        sc = mod(s-1,7) + 1;
    else
        sl = s;
        sc = s;
    end
    
    if useMarkers == 1
        option=[char(lines(sl))  colors(sc) markers(sl)]; % markers(sm) 2
    else
        option=[char(lines(sl))  colors(sc)];
    end
    if hasLineWidth == 1
        lineWidth = lineWidths(sl);
    else
        lineWidth = 2;
    end
    %option = [char(lines(sl)) colors(sc)]; 

%     i1=find(xs==1,1,'last');
%     if isempty(i1)
%      i1=1;
%     end
    i1 = 1;
    xs=xs(i1:end); 
    ys=ys(i1:end);
    ymin=max(0,min(ymin, ys(1)-0.05));
    if useMarkers == 1
        plot(xs,ys,option,'LineWidth',lineWidth,'MarkerSize',markerSize); %2
    else
        plot(xs,ys,option,'LineWidth',lineWidth); %2
    end
end

for s = 1: ns
    [xs,ys] = stairs(r(:,s),[1:np]/np);
    %i1=find(xs==1,1,'last');
    %if isempty(i1), i1=1; end
    i1 = 1;
    xs=xs(i1:end); 
    ys=ys(i1:end);
    if hastypes == 0
        sl = mod(s-1,3) + 1; 
        sc = mod(s-1,7) + 1;
    else
        sl = s;
        sc = s;
    end
    if useMarkers == 1
        options0=[char(lines(sl))  colors(sc) markers(sl)]; % markers(sm) 2
    else
        options0=[char(lines(sl))  colors(sc) markers(2)];
    end
    plot(xs(1),ys(1),options0,'MarkerFaceColor',colors(sc),'MarkerSize',6);
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

%axis([ 1 1.05*max_ratio 0 1 ]);

% Legends and title should be added.
% if ns > 3
%     %legend(legstr,'Location','SouthEast','Fontsize',7);
%     %legend(legstr,'Location','SouthEast','Fontsize',7,'interpreter','latex');
%     
% %elseif ns == 3
%     legend(legstr,'Location',leglocation,'Fontsize',8,'interpreter','latex'); % SouthEast % 10
% else
%     legend(legstr,'Location',leglocation,'Fontsize',14,'interpreter','latex');
% end

% Minimum overall ratio
xmin=min(min(r));
ymin=max(0,ymin);

if hasLimits == true
    
    xmin        = XLim(1);
    
    axis([xmin 1.05*max_ratio ymin ymax]);
    cax = gca;
    cax.XTick   = XTick;
    cax.XScale  = 'log'; 
    
    line([1 1],1.1.*[0 1],'LineStyle','--','LineWidth',0.75,'Color',[0.5,0.5,0.5])
    
    box(cax,'on');

else

    axis([xmin 1.05*max_ratio ymin ymax]);
    
end

% Legends and title should be added.
if ns > 3
    %legend(legstr,'Location','SouthEast','Fontsize',7);
    %legend(legstr,'Location','SouthEast','Fontsize',7,'interpreter','latex');
    
%elseif ns == 3
    legend(legstr,'Location',leglocation,'Fontsize',legFontSize,'interpreter','latex'); % SouthEast % 10
else
    legend(legstr,'Location',leglocation,'Fontsize',legFontSize,'interpreter','latex');
end
    
% if hasTauMax == true
% 
%     axis([xmin 1.05*tauMax ymin ymax]);
%     
% else
%     
%     axis([xmin 1.05*max_ratio ymin ymax]);
%     
% end

xlabel('$\tau$','FontSize',14,'interpreter','latex');
ylabel('$\rho_s(\tau)$','FontSize',14,'interpreter','latex');
%xlabel('\tau','FontSize',14);
%ylabel('\rho_s(\tau)','FontSize',14);

