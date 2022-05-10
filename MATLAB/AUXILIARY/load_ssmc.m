%------------------------ load_ssmc --------------------------------------%
%
% Loading the "SuiteSparse Matrix Collection"
%
% This can be used to "preload" the problems prior to running
% EXPERIMENT_II.m
%
% Note that the matrices for EXPERIMENT_I.m and EXPERIMENT_IV.m
% are already included
%-------------------------------------------------------------------------%
% 04/09/20, J.B., Initial implementation
% 05/07/22, J.B., Preparation for release

addpath '../external/ssget';

% Get index (https://sparse.tamu.edu/interfaces)

index = ssget;

% Dimensions for EXPERIMENT_II
nl = 10000;
nu = Inf;
nprob = 51;

% Dimensions for EXPERIMENT_I
% This data is already preloaded
% nl = 1000;
% nu = 10000;
% nprob   = 42;

condit =  (index.nrows > index.ncols) & ...
          (index.ncols <= nu) & ...
          (nl <= index.ncols);

% Linear indices
ids  = find(condit);
nids = length(ids);


% This can be used to preload data
for i = 1:nprob
    Prob = ssget(ids(i));
end
    
    