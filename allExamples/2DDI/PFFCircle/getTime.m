function [ current_time_in_sec ] = getTime()
%#codegen


coder.updateBuildInfo('addSourceFiles','matlab_c_get_current_time_in_sec.c');
coder.cinclude("matlab_c_timing_utils.h");

current_time_in_sec = double( 0.0 );
current_time_in_sec = coder.ceval( 'matlab_c_get_current_time_in_sec' );


end

