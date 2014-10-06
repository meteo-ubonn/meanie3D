#ifndef M3D_NAMESPACE_DEFINITIONS_H
#define M3D_NAMESPACE_DEFINITIONS_H

#include <iostream>             // namespace std
#include <netcdf>               // namespace netCDF
#include <boost/assert.hpp>     // namespace boost;

// This file manages the imports of namespaces into the 
// various namespaces defined in this software.
// That way all dependencies can be changed in one place
// should the need arise. Also, this is a great way of 
// actually seeing all namespaces and dependencies.

extern "C"
{
    namespace m3D 
    {
        namespace utils {
            namespace netcdf {}
            namespace sets {}
            namespace vectors {}
        }
    }
}

#endif
