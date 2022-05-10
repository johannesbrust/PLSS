function x =iter_GaussBpd(A,b,x,options)
eta =randn(length(b),options.block_size); %normrnd(0,1, [size(b),options.block_size] );
Aeta = A*eta;
etaAeta = eta'*Aeta;
x = x - eta* (etaAeta\(Aeta'*x -(eta')*b));
end