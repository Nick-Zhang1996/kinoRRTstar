% clc;
% close all;
% clear all;

function [retval] = benchmarkRRT()

%tic
num_of_runs = 1;
run_RRTstar = 1;

dim=2;
 
if dim == 2
    stepsize = 4;
    radius = 3;
    samples = 4000;
elseif dim == 3
    stepsize = 3;
    radius = 3;
    samples = 1000;
end

random_world = 0;
show_output = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sits = 1:size(stepsize,2)
    
    segmentLength = stepsize(sits);

    if run_RRTstar == 1
    avg_its = 0;
    avg_path = 0;
        for i = 1:num_of_runs
            [n_its, path_n] = RRTstar3D(dim, segmentLength, radius, random_world, show_output, samples);
            avg_its = avg_its + n_its;
            avg_path = avg_path + path_n;
        end

    end
    
end
retval = 0;
%toc
end