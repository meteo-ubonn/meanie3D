#ifndef M3D_TIME_UTILS_H
#define M3D_TIME_UTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <sys/time.h>
#include <iostream>

namespace m3D { namespace utils { 

    timeval start_time_;
    void start_timer(const std::string& message = "")
    {
        if (!message.empty())
        {
            cout << message << endl;
            fflush(stdout);
        }
        gettimeofday(&start_time_,NULL);
    }

    double stop_timer()
    {
        timeval end_time;
        gettimeofday(&end_time,NULL);

        return double(end_time.tv_sec-start_time_.tv_sec)+ double(end_time.tv_usec-start_time_.tv_usec)/1000000;
    }
}}

#endif
