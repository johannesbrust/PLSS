%setup_RandomLinearLab
%
%
%    Adds to path the files of the current folder

% Get root location of resultsSetup
root = fileparts(which(mfilename)); 

% ----------------------------------------------------------------------
% Add the appropriate subdirs to path.
% ----------------------------------------------------------------------
addpath(genpath(root));
fprintf('%s \n',root);
fprintf(['\n The above directories and their subdirectories have been'  ,...
         ' temporarily added to your path.\n']);     
     
clear;
