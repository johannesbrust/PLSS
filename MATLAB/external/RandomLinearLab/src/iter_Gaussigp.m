function x =iter_Gaussigp(A,b,x,options)
eta = randn(size(b));
ATeta = A'*eta;
x = x -eta'*( A*x -b)*(ATeta./(norm(ATeta)^2));
end