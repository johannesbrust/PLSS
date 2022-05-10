function x= proxy_sample(A,b,x,invB,S)
AtS = A'*S; 
invBAtS = invB *(AtS);
x = x - invBAtS*inv((AtS')*invBAtS)*(AtS'*x-S'*b);
end