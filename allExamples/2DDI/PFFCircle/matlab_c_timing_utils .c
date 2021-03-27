#if defined( _MATLAB_C_TIMING_UTILS_C_ )
    /* This file has already been included. */
#else

#  define _MATLAB_C_TIMING_UTILS_C_

#  include "matlab_c_timing_utils.h"

/*
 * matlab_c_get_current_time_in_sec: returns the current time (that is, the
 * time since the application started) in seconds. In order to only return
 * MATLAB-friendly data types, this function converts the time to seconds
 * and produces a "double". However, in some cases, the results may be
 * smaller than "double"'s epsilon.
 */
double matlab_c_get_current_time_in_sec( void )
{

    const clock_t current_time_in_system_units = clock();
    
    const double current_time_in_s = 
            ( (double) current_time_in_system_units )/CLOCKS_PER_SEC;

    return( current_time_in_s );

}

#endif /*  _MATLAB_C_TIMING_UTILS_C_ */
