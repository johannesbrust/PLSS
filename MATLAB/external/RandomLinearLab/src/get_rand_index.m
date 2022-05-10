function [j,options] = get_rand_index(options)
% Test if random block has been used up
if(options.rand_index> options.rand_block_size)
    options.S= gendist(options.pna,options.rand_block_size,1);
    options.rand_index =1;
end
% Pick up a random index
j = options.S(options.rand_index);
options.rand_index = options.rand_index+1;
end