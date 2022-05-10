function x =iter_LSCD(A,b,x,options)
% Update Ax -b, don't recalculate.
[j,options] = get_rand_index(options);
% calculate descent direction
d  =(A(:,j)'*(options.grad)/options.na(j));
x(j) = x(j) -d;
% update gradient
options.grad = options.grad -A(:,j)*d;
assignin('caller', 'options', options);
end