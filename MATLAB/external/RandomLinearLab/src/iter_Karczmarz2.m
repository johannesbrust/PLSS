function x =iter_Karczmarz2(A,b,x,options)
%[j,options] = get_rand_index(options);
j = randi(size(A,1));
x = x -(((A(j,:)*x -b(j))/norm(A(j,:))^2)).*A(j,:)';
assignin('caller', 'options', options);
end
