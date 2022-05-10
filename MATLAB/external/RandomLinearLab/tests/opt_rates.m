%% selects the optimal probabilities for given sampling vectors [s_1,\ldots, s_n] = S.
%% Peter Richtarik  - 2015
function [p, Value_OPT,optp_time] = opt_rates(A,B,S)
tic;

V = inv(full(B)^(1/2))*A'*S;
% Problem: max_{p} [ t  = lambda_min (sum_i p_i v_i*v_i^T) ],
% where v_i are unit norm vectors and p is a probability vector

% We know that optimal t satisfies these bounds:
% 0 <= t <= 1/n


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a random problem, so not very good for testing
sA = size(A);
m = sA(1);
n = sA(2);
r= size(S,2);
%V = rand(n,n); % V = [v_1, ..., v_n]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE OPTIMAL PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(r,1); 
s = diag(V'*V);
D = diag(s);
W = V*inv(D^(1/2)); % columns of W are now of unit norm

cvx_quiet false

cvx_begin sdp
    variables t p(r,1) 
    maximize t
    subject to
        S = W*diag(p)*W';
        e'*p == 1;
        p >= 0;
        S >= t*eye(n); % we want the minimal eigenvalue of S to be at least t
cvx_end
optp_time = toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVENIENT PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = s/sum(s); % convenient probabilties
Value_CONVENIENT = min(eig(W*diag(q)*W'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Value_OPT = cvx_optval;

display(sprintf('Convenient eigenvalue = %2.3e',Value_CONVENIENT))
display(sprintf('Optimal eigenvalue    = %2.3e',Value_OPT))
display(sprintf('Upper bound ( = 1/n)  = %2.3e', 1/n))
display(sprintf('Time  = %2.3f', optp_time))
end
