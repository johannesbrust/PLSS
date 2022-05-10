function output_stats = distribution_output(Prob, iter_func, boot_func,options,num_trials, pq )

TIMES =[];
ERRORS =[];
%% Performing multiple tests, collecting data
for i =1:num_trials
    [~, output]= solve_system(Prob,iter_func, boot_func,options );
    lent(i) = length(output.times);
    %set the max iterations for uniform plot sizez
    if(i==1)
        options.max_iterations = lent(i)-1;
    end
    
    try
        TIMES = [TIMES ; output.times] ;
        ERRORS = [ERRORS; output.errors] ;
    catch errooerer
        errooerer
        pause;
    end
    
    if(i==1) % Must hit max iterations
        options.max_time =10^6;
        options.ep=0;
    end
end

%% Thining 
        lent = length(TIMES);
        if(lent>2000)
            grid =1:floor(lent/2000):lent;
          %  lent=length(grid);
            TIMES = TIMES(:,grid);
            ERRORS = ERRORS(:,grid);
        end
%% Calculatings stats

output_stats = output;

    TIMES_quant = zeros(2,length(TIMES));
    TIMES_mean = zeros(1,length(TIMES));
    ERRORS_quant =zeros(2,length(TIMES));
    ERRORS_mean = zeros(1,length(TIMES));
    for i=1:length(TIMES)
        TIMES_quant(:,i) = [quantile(TIMES(:,i),1-pq );  quantile(TIMES(:,i),pq)];
        ERRORS_quant(:,i) = [quantile(ERRORS(:,i),1-pq) ;  quantile(ERRORS(:,i),pq)];
        TIMES_mean(i) =  quantile(TIMES(:,i),0.5);
        ERRORS_mean(i) = quantile(ERRORS(:,i),0.5);
    end
    output_stats.times = TIMES_mean;
    output_stats.errors =ERRORS_mean;
    output_stats.time_stats= TIMES_quant;
    output_stats.error_stats= ERRORS_quant;
  %  TIMES_quant =TIMES_quant;
  %  ERRORS_quant =ERRORS_quant;
        
end