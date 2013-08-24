#ifndef _M3D_Matrix_H_
#define _M3D_Matrix_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

namespace m3D { namespace utils {

    // Matrix
    
    template <typename T>
    struct Matrix
    {
        // A matrix is a 2D vector construct in this context
        typedef vector< vector<T> > matrix_t;
        
        /** Creates a new matrix of the given dimensions
         */
        static matrix_t
        create_matrix(size_t width, size_t height);
        
        // Flag Matrix
        
        typedef vector< vector<int> > flag_matrix_t;
        
        /** Creates a new matrix of the given dimensions
         */
        static flag_matrix_t
        create_flag_matrix(size_t width, size_t height, int defaultValue=0 );
    };
    
}}

#endif
