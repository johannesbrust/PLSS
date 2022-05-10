function Prob = uniform_Prob(m,n)
Prob.A = randn(m,n);   % randn(n,n)
Prob.sol = randn(n,1);
Prob.b = Prob.A*Prob.sol;
Prob.title =[ 'uniform-random' num2str(m) 'X' num2str(n)];
end