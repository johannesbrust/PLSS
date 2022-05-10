function x =iter_Gausspd(A,b,x,options)
eta = randn(size(b));
Aeta = A*eta;
etaAeta = eta'*Aeta;
x = x -( Aeta'*x -(eta')*b)*(eta./etaAeta);
end