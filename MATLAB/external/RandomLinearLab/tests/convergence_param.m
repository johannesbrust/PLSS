function cnv_param= convergence_param(Prob,method_name,method_number)

A = Prob.A;
if(~isempty(method_name))
    if (issparse(A))
        Asvd=svds(A);
    else
        Asvd=svd(A);
    end
    switch method_name      
        case 'Gauss pd'
            cnv_param = min(Asvd)/sum(Asvd);
            % Based on Conjecture
            %cnv_param =sqrt(min(eig(A)))/trace(sqrt(A));
        case 'Gauss LS'
            cnv_param =  min(Asvd)^2/ (norm(A,'fro')^2);
            %  cnv_param =  min(svd(A))^2/ (length(Prob.sol)*max(svd(A))^4);
            % Based on Conject
            %Asvd=svd(A);
            %cnv_param = min(Asvd)/sum(Asvd);
        case 'Gauss Kaczmarz'
            cnv_param =  min(Asvd)^2/ (norm(A,'fro')^2);
            % Based on Conject
            %Asvd=svd(A);
            %cnv_param = min(Asvd)/sum(Asvd);
        case 'CD pd'
            cnv_param =min(Asvd)/trace(A);
        case 'CD LS'
            cnv_param =  min(Asvd)^2/ (norm(A,'fro')^2);
        case 'Kaczmarz'
            cnv_param = min(Asvd)^2/ (norm(A,'fro')^2);
        otherwise
            display('convergence_param: I do not know that method_name');
            cnv_param =-2;
    end
end
end
%(log(norm(options.sol)^2/options.ep))); 
% Gaussls
%  ceil(trace(sqrt(A))/sqrt(min(eig(A))) 
% 
% Gausspd
% ceil( (min(Asvd)/sum(Asvd)))  
%  
% Gaussigp 
%  ceil( (min(Asvd)/sum(Asvd))) 
%  
%  
% CD
%  ceil(trace(A)/min(eig(A)))
%  
% LSCD
% 
% ceil((norm(A,'fro')^2/min(eig(A))^2) 
% 
% 
% Karczmarz
% 
% ceil((norm(A,'fro')^2/min(eig(A))^2)