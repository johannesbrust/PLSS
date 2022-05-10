function [Prob] = sparse_matlab(m,n,density,rc)
Prob.A = sprand(m,n,density,rc) ;  % sprandsym(n,density,rc,kind)
Prob.title =['sprandn-' num2str(m) '-' num2str(n)  '-' num2str(density) '-' num2str(rc)  ];
Prob.sol = rand(n,1);
Prob.b = Prob.A*Prob.sol;
% Prob.b = rand(m,1);
% Prob.sol = Prob.A\Prob.b;
end