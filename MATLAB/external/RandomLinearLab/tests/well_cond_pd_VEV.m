function Prob =well_cond_pd_VEV(n)
V = randn(n,n); 
[Q,~] = qr(V);
eigsn = abs(rand(n,1));
Prob.A= Q*diag(eigsn)*Q';
Prob.b = rand(n,1);
Prob.sol = Prob.A\Prob.b;
Prob.title =['random posdef ' num2str(n) 'X' num2str(n)  ' with cond= ' num2str(cond(Prob.A)) ];
end