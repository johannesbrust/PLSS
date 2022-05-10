function output_theos = theoretical_output(Prob,output_reals)

output_theos ={};
if(~iscell(output_reals))
    output_reals = {output_reals};
end
for j = 1:length(output_reals)
    output_real = output_reals{j};
    lt = length(output_real.times);
    output_theo.times = 1:1:lt;
    output_theo.errors(1) =100;
    C= convergence_param(Prob,output_real.name,[]);
    for i =2:lt
        output_theo.errors(i) = (1-C)*output_theo.errors(i-1);
    end
    output_theo.name =['theo. ' output_real.name ];
    output_theos = [output_theos; output_theo];
end

end

