#ifndef M3D_RANDOM_UTILS_H
#define M3D_RANDOM_UTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

namespace m3D { namespace utils { 

    /** @returns uniform random variable (0..1)
     */
    float ranf();

    /** Normal random variate generator, courtesy of
     * ftp://ftp.taygeta.com/pub/c/boxmuller.c
     * @param m mean value
     * @param s standard deviation
     * @return random variable with mean m, standard deviation s
     */
    float box_muller(float m, float s);
}}

#endif
