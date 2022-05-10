%% fix PRob 

if( tril( full(Prob.A) ) ==full(Prob.A) )
    Prob.A = Prob.A +Prob.A' -diag(diag(Prob.A));
    spy(Prob.A); 
end
if(~isfield(Prob,'sol'))
    SA =size(Prob.A);
    Prob.sol = rand(SA(1),1);
    Prob.b = Prob.A*Prob.sol;
end