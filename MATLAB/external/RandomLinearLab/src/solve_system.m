% A wrapper function for testing and timing iterative methods for
% solving a linear system. Note that only the iteration function 'iter_func'
% and the boot function 'boot_func' are timed - 2015 - Robert M. Gower
%--------------------------------------------------------------------------
% 07/02/20, J.B., Modification to include errorScale option for comparisons
function [x, output] = solve_system(Prob, iter_func, boot_func,options )
A = Prob.A;
b= Prob.b;
tic;
x = boot_func(A,b,options);
times(1) = toc;
if isfield(options,'errorScale') 
    errorScale = options.errorScale;
else
    errorScale = 100; 
end
if isfield(options,'errorSqrt') 
    errorSqrt = options.errorSqrt;
else
    errorSqrt = 0; 
end
err = calculate_error(x,Prob,options);
errors(1) = errorScale*err; % Initial relative error
if errorSqrt == 1
    errors(1) = sqrt(errors(1));
end
for i = 1:options.max_iterations 
%{
    if (mod(i,30)==1) %prints a header ever 30 iterations
        fprintf('-------------------\n'); fprintf('Iter | Error   |  Time   \n'); fprintf('-------------------\n');
    end
    fprintf('%d  | %3.4f  |  %3.4f \n',i,errors(i),times(i) );
%}
    if(times(i) >options.max_time )
        break;
    end
    tic;
    x = iter_func(A,b,x,options);  % Computes x^(k+1) from x^k, A and b.
    times(i+1)= times(i) +  toc;
    err = calculate_error(x,Prob,options);
    errors(i+1) = errorScale*err; % Calculate relative error
    if errorSqrt == 1
        errors(i+1) = sqrt(errors(i+1));
        err = errors(i+1);
    end
    if(err < options.ep)
        fprintf('%d  | %3.4f  |  %3.4f \n',i,errors(i+1),times(i+1) );
        break;
    end
end
output.times = times;
output.errors = errors;
output.name = options.name;
end
