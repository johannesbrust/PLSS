function [Prob] = sparse_symmetric_matlab(n,density,rc,type)

if(isempty(type))
    Prob.A = sprandsym(n,density,rc) ;  % sprandsym(n,density,rc,kind)
    Prob.title =['sprandsym-' num2str(n)  '-' num2str(density) '-' num2str(rc)  ];
else
    Prob.A = sprandsym(n,density,rc,type);
    Prob.title =['sprandsym-' num2str(n) '-' num2str(density) '-' num2str(rc) '-' num2str(type)  ];
end
Prob.b = rand(n,1);
Prob.sol = Prob.A\Prob.b;
end