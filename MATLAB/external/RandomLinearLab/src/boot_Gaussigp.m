function x = boot_Gaussigp(A,b,options)
sA = size(A);
options.iter=(sA(1)*sA(2))^2;
% Theoretical number of iterations required to reach ep accuracy in
% expectation
% Asvd= svd(A);
%options.iter =  ceil( (min(Asvd)/sum(Asvd)) *(log(norm(options.sol)^2/options.ep))); 
options.name= 'Gauss Kaczmarz';
if(~isfield(options,'metric'))
   options.metric = eye(sA(2));
end
assignin('caller', 'options', options);
x =zeros(sA(2),1);
end