#ifndef M3D_OPERATION_ITERATE_H
#define M3D_OPERATION_ITERATE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace.h>

namespace m3D { 

    template <typename T>
    class IterationOperation : public Operation<T>
    {

    public:

        IterationOperation( FeatureSpace<T> *fs, PointIndex<T> *index ) : Operation<T>( fs, index ) {};

        virtual ~IterationOperation() {}

        /** Performs a mean-shift iteration from the given starting
         * point until one of the two termination criteria is met.
         * 
         * @param origin of the iteration (a point in feature-space)
         * @param params 
         * @param kernel
         * @param weight weight function to use (or NULL)
         * @param termcrit_epsilon smallest iteration step
         * @param termcrit_iterations maximum number of iterations
         */
        void
        iterate(Point<T> *origin, 
                const SearchParameters *params, 
                const Kernel<T> *kernel, 
                const WeightFunction<T> *weight,
                const T termcrit_epsilon,
                const size_t termcrit_iterations);

        /** Identical to 'iterate' with the difference that the iteration
         * trajectory is recorded and returned.
         * 
         * @param origin of the iteration (a point in feature-space)
         * @param params 
         * @param kernel
         * @param weight weight function to use (or NULL)
         * @param termcrit_epsilon smallest iteration step
         * @param termcrit_iterations maximum number of iterations
         * 
         * @return The trajectory of the iteration
         */
        typename FeatureSpace<T>::Trajectory * 
        get_trajectory(Point<T> *origin, 
                       const SearchParameters *params, 
                       const Kernel<T> *kernel, 
                       const WeightFunction<T> *weight,
                       const T termcrit_epsilon,
                       const size_t termcrit_iterations);


        
    };   
}

#endif
