function x =iter_BCD(A,b,x,options)
% 04/02/20, J.B., Modification to use the function randperm
% instead of randsample, b/c randperm is included with the 
% standard Matlab distribution, while randsample uses the
% Statistics and Machine Learning Toolbox.
s = randperm(length(b),options.block_size);
%s = randsample(length(b),options.block_size);
x(s) = x(s) -(A(s,s)\(A(s,:)*x -b(s))) ;
end


