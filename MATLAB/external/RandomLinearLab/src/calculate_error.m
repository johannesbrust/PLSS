function error =  calculate_error(x,Prob,options)

metric_norm = ((x-Prob.sol)'*options.metric)*(x-Prob.sol);
metric_initial = ((Prob.sol)'*options.metric)*(Prob.sol);
if(metric_initial ==0)
    error = metric_norm;
    return;
end
error = metric_norm/metric_initial;

end