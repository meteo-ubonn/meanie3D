#ifndef M3D_UTILS_TIMER_H
#define M3D_UTILS_TIMER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <sys/time.h>
#include <iostream>

namespace m3D { namespace utils { 

    timeval start_time_;

    /** Starts the timer, optionally with a message.
     * 
     * @param message
     */    
    void start_timer(const std::string& message = "")
    {
        if (!message.empty())
        {
            cout << message << endl;
            fflush(stdout);
        }
        gettimeofday(&start_time_,NULL);
    }

    /** Stops the timer and returns the seconds since
     * the last call of start_timer
     * 
     * @return 
     */
    double stop_timer()
    {
        timeval end_time;
        gettimeofday(&end_time,NULL);

        return double(end_time.tv_sec-start_time_.tv_sec) 
                + double(end_time.tv_usec-start_time_.tv_usec) / 1000000.0;
    }
}}

#endif
