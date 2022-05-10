function x = boot_Karczmarz(A,b,options)
sA = size(A);
[options.pna,options.na]= row_weight_prob(A);
options=boot_random_weighted_selection(options,sA(1));
options.name= 'Kaczmarz';
if(~isfield(options,'metric'))
   options.metric = eye(sA(2));
end
% Get probabilities
if(isfield(options,'probs')) % User supplied
    options.pna =options.probs';
    if(isequal(ones(length(options.pna))*options.pna(1),options.pna))
        options.name = [options.name '-puniforme'];
    else
        options.name =  [options.name '-popt'];
    end
end
% Theoretical number of iterations needed
%options.iter =  ceil((norm(A,'fro')^2/min(eig(A))^2)*(log(norm(options.sol)^2/options.ep))); 
assignin('caller', 'options', options);
x =zeros(sA(2),1);
end
