#ifndef M3D_OPERATION_ITERATE_IMPL_H
#define M3D_OPERATION_ITERATE_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include "iterate_op.h"

namespace m3D {

#pragma mark -
#pragma mark Meanshift Iteration    

    template <typename T> 
    typename FeatureSpace<T>::Trajectory * 
    IterationOperation<T>::iterate( Point<T> *origin, 
                                    const SearchParameters *params,                                     
                                    const Kernel<T> *kernel, 
                                    const int weight_index,
                                    const T termcrit_epsilon,
                                    const size_t termcrit_iterations )
    {
        // Allocate a fresh trajectory

        typename FeatureSpace<T>::Trajectory *trajectory = new typename FeatureSpace<T>::Trajectory();

        // Add the origin as first point

        vector<T> x = origin->values;

        trajectory->push_back( x );

        // Go

        T dx = std::numeric_limits<T>::max();

        vector<T> d( x.size() );

        size_t iter = 0;

        while ( iter < termcrit_iterations && dx >= termcrit_epsilon ) 
        {
#if DEBUG_ITERATION
            std::cout << "Iteration " << iter << " from " << x;
#endif
            // get the mean-shift

            MeanshiftOperation<T> op( this->feature_space, this->point_index );

            vector<T> shift = op.meanshift( x, params, kernel, weight_index );

            if ( iter == 0 )
            {
                origin->shift = shift;
            }

#if DEBUG_ITERATION
            std::cout << " with mean-shift " << shift;
#endif            

            // calculate termination criteria
            dx = (T)vector_norm( shift );

            // calculate iteration end point (re-use shift variable)
            // Use a slightly optimized loop form 

            size_t index = 0; 
            typename vector<T>::iterator it;
            for ( it=shift.begin(); it!= shift.end(); it++ )
            {
                *it += x[index++];
            }

#if DEBUG_ITERATION
            std::cout << " to " << shift;
#endif            
            // add to trajectory

            trajectory->push_back( shift );

#if DEBUG_ITERATION
            std::cout << " (dx = " << dx << ")" << std::endl;
#endif            
            iter++;

            x = shift;
        }

        // substract one, since termination criterion is checked
        // at beginning of the loop
        iter--;

        return trajectory;
    }
}

#endif
