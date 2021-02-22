cfg = coder.config('exe')
codegen -report -config cfg benchmarkRRT main.c main.h
system("benchmarkRRT.exe")
% load('H:\DCSL\dongliang_rrt_star_v2\codegen\lib\benchmarkRRT\buildInfo.mat')
% packNGo(buildInfo, 'fileName', 'benchmarkRRT.zip');