function x =iter_CD(A,b,x,options)
[j,options] = get_rand_index(options);
x(j) = x(j) -((A(j,:)*x -b(j))/options.na(j));
assignin('caller', 'options', options);
end