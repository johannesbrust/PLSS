%% Options parameters
OUTPUTS ={};
options.ep=10^(-3);  % impossible  precison, must run untill time depleets. 
options.max_time = 0.5;  % 10 min
options.max_iterations =  min(sum(size(Prob.A))^5, 10^10);
options.metric = Prob.A'*Prob.A;   % this forces all methods to use the same error measurement
%% The quantile percentage used
pq = 0.05;
nm =4;
lent = zeros(1,nm);
TIMES =cell(1,nm);
ERRORS=cell(1,nm);
num_trials = 10;
max_iterations = ones(1,nm)*options.max_iterations;
%% Performing multiple tests, collecting data
for i =1:num_trials
    OUTPUTS ={};
    OUTPUTS= all_overdet_methods(OUTPUTS, options, Prob,nm,max_iterations);
    for j =1:length(OUTPUTS)
        
        %set the max iterations for uniform plot sizez
        if(i==1)
            lent(j) = length(OUTPUTS{j}.times);
            max_iterations(j) = lent(j)-1;
        end
        
        try
            TIMES{j} = [TIMES{j} ; OUTPUTS{j}.times] ;
            ERRORS{j} = [ERRORS{j}; OUTPUTS{j}.errors] ;
        catch errooerer
            errooerer
            pause;
        end
        
    end
    if(i==1) % Must hit max iterations
        options.max_time =10^6;
        options.ep=0;
    end
end
%% Thining the set
    for j =1:length(TIMES)
        lent(j) = length(TIMES{j});
        if(lent(j)>2000)
            grid =1:floor(lent(j)/2000):lent(j);
            lent(j)=length(grid);
            TIMES{j} = TIMES{j}(:,grid);
            ERRORS{j} = ERRORS{j}(:,grid);
        end
    end
%%
TIMES_quant =cell(1,nm);
ERRORS_quant=cell(1,nm);
line_colors= lines(nm);
OUTPUTS_mean= OUTPUTS;
for i =1:length(TIMES)
    timesi = TIMES{i};
    timesi_quant =[];
    timesi_mean =[];
    errorsi = ERRORS{i};
    errorsi_quant =[];
    errorsi_mean = [];
    for j=1:length(timesi)
        timesi_quant(:,j) = [quantile(timesi(:,j),1-pq );  quantile(timesi(:,j),pq)];
        errorsi_quant(:,j) = [quantile(errorsi(:,j),1-pq) ;  quantile(errorsi(:,j),pq)];
        timesi_mean(j) =  quantile(timesi(:,j),0.5);
        errorsi_mean(j) = quantile(errorsi(:,j),0.5);
    end
    OUTPUTS_mean{i}.times = timesi_mean;
    OUTPUTS_mean{i}.errors =errorsi_mean;
    TIMES_quant{i} =timesi_quant;
    ERRORS_quant{i} =errorsi_quant;
    X = [timesi_mean,fliplr(timesi_mean)];
    Y = [errorsi_quant(1,:) , fliplr(errorsi_quant(2,:))]; palpha=0.4;
    h=fill(X,Y,line_colors(i,:),'FaceAlpha',palpha,'EdgeAlpha',palpha);
    set (gca, 'Yscale', 'log');
    set(h,'EdgeColor','None');
    set(h,'facealpha',.2);
    %alpha(.4)
    hold on;
end
Probplot.title = [Prob.title '-dist'];
prettyPlot_solve_system_wrapper(OUTPUTS_mean,Probplot)

