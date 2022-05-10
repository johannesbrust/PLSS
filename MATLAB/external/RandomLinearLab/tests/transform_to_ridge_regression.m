function Prob= transform_to_ridge_regression(Prob)
sA = size(Prob.A);
if(~isempty(findstr(Prob.title, '-ridge')))
    return;
end
% Making a least-squares system if it's not square
if(sA(1) ~= sA(2))
   if(isfield(Prob,'b'))
    Prob.b= (Prob.A)'*Prob.b;
   end
   Prob.A = (Prob.A)'*Prob.A; 
end

% Adding regularization paramenter
Prob.A = Prob.A+eye(sA(2));
Prob.title=[Prob.title '-ridge']
end