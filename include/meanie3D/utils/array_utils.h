#ifndef M3D_UTILS_DEBUG_H
#define M3D_UTILS_DEBUG_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <stdlib.h>
#include <iostream>

namespace m3D { namespace utils { 

    /** Prints the given array.
     */
    template <typename T> 
    void print_array( T *array, size_t count )
    {
        using namespace std;
        
        cout << "(";

        for ( size_t i = 0; i < count; i++ )
        {
            cout << array[i];

            if ( i < count - 1 )
            {
                cout << ",";
            }
        }       

        cout << ")";
    }
}}

#endif


