
% clc;
% close all;
% clear all;

function [retval] = benchmarkRRT()
coder.extrinsic('tic');
coder.extrinsic('toc');
tic
num_of_runs = 1;
run_RRTstar = 1;

dim=3;
 
if dim == 2
    stepsize = 4;
    radius = 3;
    samples = 2000;
elseif dim == 3
    stepsize = 4;
    radius = 4;
    samples = 600;
end

show_output = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    segmentLength = stepsize;

    if run_RRTstar == 1

        for i = 1:num_of_runs
            RRTstar3D(dim, segmentLength, radius, show_output, samples);
        end

    end
    
retval = 0;
toc
end
