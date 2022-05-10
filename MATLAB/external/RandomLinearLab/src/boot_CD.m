function x = boot_CD(A,b,options)
sA = size(A);
options.name= 'CD pd';
options.grad = -b;
if(~isfield(options,'metric'))
    options.metric = A;
end
% Get norm diagonal, needed for iterates.
for i = 1: sA(1)
    options.na(i) = norm(A(i,i));
end
% Get probabilities
if(isfield(options,'probs')) % User supplied
    options.pna =options.probs';
    if(isequal(ones(length(options.pna))*options.pna(1),options.pna))
        options.name = [options.name '-puniforme'];
    else
        options.name =  [options.name '-opt'];
    end
else  % Convenient probabilities proportional to diagonal
    options.pna= options.na/sum(options.na);
end
options=boot_random_weighted_selection(options,sA(1));
%options.iter =  ceil(trace(A)/min(eig(A))) *(log(norm(options.sol)^2/options.ep)); % Theoretical number of iterations

assignin('caller', 'options', options);
x =zeros(sA(2),1);
end