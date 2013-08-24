#ifndef _VERBOSITY_H_
#define _VERBOSITY_H_

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