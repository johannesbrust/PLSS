function x = boot_LSCD(A,b,options)
sA = size(A);
[options.pna,options.na]= row_weight_prob(A'); % column weights because of A'
options=boot_random_weighted_selection(options,sA(1));
%options.iter =  ceil((norm(A,'fro')^2/min(eig(A))^2) *(log(norm(options.sol)^2/options.ep)));
options.name=  'CD LS';
x =zeros(sA(2),1);
options.grad = -b;
if(~isfield(options,'metric'))
   options.metric = A'*A;
end
assignin('caller', 'options', options);

end