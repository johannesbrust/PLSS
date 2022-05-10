function x = boot_BCD(A,b,options)
sA = size(A);
if(~isfield(options,'block_size'))
    options.block_size =ceil(sqrt(sA(2)));
end
options.name= 'Block CD-pd';
if(~isfield(options,'metric'))
   options.metric = A;
end
assignin('caller', 'options', options);
x =zeros(sA(2),1);
end