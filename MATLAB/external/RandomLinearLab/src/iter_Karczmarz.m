function x =iter_Karczmarz(A,b,x,options)
[j,options] = get_rand_index(options);
x = x -(((A(j,:)*x -b(j))/options.na(j))).*A(j,:)';
assignin('caller', 'options', options);
end