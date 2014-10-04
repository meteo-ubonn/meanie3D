#ifndef M3D_VERBOSITY_H
#define M3D_VERBOSITY_H

namespace m3D { namespace utils {

    /** Verbosity
     */
    typedef enum {
        VerbositySilent = 0,
        VerbosityNormal = 1,
        VerbosityDetails = 2,
        VerbosityAll = 3
    } Verbosity;

}}

#endif