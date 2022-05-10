function options=boot_random_weighted_selection(options,rand_block_size)
options.rand_block_size = rand_block_size;
options.rand_index =1;
options.S= gendist(options.pna,options.rand_block_size,1);
end