function Prob = Gauss_pd_Prob(n)
A = randn(n,n);   % randn(n,n)
Prob.A = A'*A;
Prob.b = rand(n,1);
Prob.sol = Prob.A\Prob.b;
Prob.title =[ 'uniform-random-' num2str(n) 'X' num2str(n)];
end