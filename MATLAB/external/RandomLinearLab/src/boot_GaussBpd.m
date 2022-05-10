function x = boot_GaussBpd(A,b,options)
sA = size(A);
options.iter=(sA(1)*sA(2))^2;
if(~isfield(options,'block_size'))
    options.block_size =ceil(sqrt(sA(2)));
end
options.name= 'Block Gaussain pd ';
if(~isfield(options,'metric'))
   options.metric = A;
end
assignin('caller', 'options', options);
x =zeros(sA(2),1);
end