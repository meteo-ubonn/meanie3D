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

        /** Iterative mean-shift, starting at position x in feature-space. 
         * @param x iteration origin
         * @param bandwidth vector representing the diagonal of the bandwidth matrix
         * @param kernel 
         * @param weighing function : externally handed in weighing function, can be NULL
         * @param termcrit_epsilon : if the distance between iterations gets lower or equal to this value,
         *                           the iteration is stopped.
         * @param termcrit_iterations : if more than this number of iterations have been done, it stops.
         * @return The trajectory of the iteration
         */
        typename FeatureSpace<T>::Trajectory * 
        iterate( Point<T> *origin, 
                 const SearchParameters *params, 
                 const Kernel<T> *kernel, 
                 const int weight_index,
                 const T termcrit_epsilon,
                 const size_t termcrit_iterations );

    };   
}

#endif
