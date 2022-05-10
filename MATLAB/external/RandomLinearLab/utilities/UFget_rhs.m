cd /home/s1065527/Software/Matlab/UFget

index = UFget;
f = get_rhs_num();
[y, j] = sort (index.nrows (f)) ;
f = f (j) ;
for i = 1:length(f)
    Prob{i} = UFget (f(i)) ;
    fprintf ('%4d %s\n', i, Prob{i}.name) ;
end

%%
index = UFget;
i_posdef=  find(index.posdef);
i_fposdef = intersect(i_posdef,f);
[y, j] = sort (index.nrows (i_fposdef)) ;
i_fposdef = i_fposdef (j) ;

count =1;
for i = 1:length(Prob)
      if (isfield (Prob{i}, 'posdef'))
         fprintf ('%4d %s\n', i, Prob{i}.name) ;
         Prob_posdef{count} =Prob{i};
         count= count+1;
     end  
end
