%% Add Path 

if ispc
    addpath('C:\gurobi651\win64\matlab\');
else
    addpath('~/software/gurobi651/linux64/matlab/');
end
gurobi_setup
addpath(genpath(pwd))
