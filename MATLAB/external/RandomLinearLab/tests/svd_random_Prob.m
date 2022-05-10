function Prob =svd_random_Prob(m,n)
V = randn(n,n); 
U = randn(m,m);
[Q,R] = qr(V);
[QU,RU] = qr(U);
eigsn = abs(rand(min(m,n),1));
Prob.A= (Q*[diag(eigsn) zeros(n,m-n)]*QU')';
Prob.sol =  rand(n,1);
Prob.b = Prob.A*Prob.sol;
%%inconsistent version
%Prob.b = rand(m,1);
%Prob.sol = Prob.A\Prob.b;
Prob.title =  ['random\_svd-' num2str(m) '-' num2str(n) ];
end