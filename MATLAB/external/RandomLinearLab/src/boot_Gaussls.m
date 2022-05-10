function x = boot_Gaussls(A,b,options)
sA = size(A);
options.iter=(sA(1)*sA(2))^2;
% Theoretical number of iterations required to reach ep accuracy in
% expectation
%options.iter =  ceil(trace(sqrt(A))/sqrt(min(eig(A))) *(log(norm(options.sol)^2/options.ep))); 
options.name= 'Gauss LS';
if(~isfield(options,'metric'))
   options.metric = A'*A;
end
assignin('caller', 'options', options);
x =zeros(sA(2),1);
end