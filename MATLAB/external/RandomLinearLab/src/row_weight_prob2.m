function [pna,na]= row_weight_prob2(A)
%na = diag(A*A');
%{
sA  =size(A);
for i = 1: sA(1)
    na(i) = norm(A(i,:))^2;
end
pna= na/sum(na);
%}
na  = ones(1,size(A,1));
pna = na;

end
