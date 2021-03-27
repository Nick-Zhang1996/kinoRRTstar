%% generate exe
cfg = coder.config('exe')
codegen -report -config cfg benchmarkRRT main.c main.h
fprintf("running program... \n")
system("benchmarkRRT.exe")

%% package
%load('H:\DCSL\dongliang_rrt_star_v2\codegen\lib\benchmarkRRT\buildInfo.mat')
%packNGo(buildInfo, 'fileName', 'D:\Documents\GitHub\kinoRRTstar\allExamples\2DDI\PFFCircle\benchmarkRRT.zip');