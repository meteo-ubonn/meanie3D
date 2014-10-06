#ifndef M3D_OPERATION_MEANSHIFTOPERATION_H
#define M3D_OPERATION_MEANSHIFTOPERATION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/operations.h>
#include <meanie3D/featurespace.h>
#include <meanie3D/index.h>

#include <vector>

namespace m3D { 

    template <typename T>
    class MeanshiftOperation : public Operation<T>
    {
    public:

        //static const int NO_WEIGHT;

        MeanshiftOperation(FeatureSpace<T> *fs,
                           PointIndex<T> *index ) : Operation<T>(fs,index) {}

        virtual ~MeanshiftOperation() {}


        /** Call this up-front to preempt lazy index construction.
         * @param search parameters
         */
        void prime_index(const SearchParameters *params);

        /** Meanshift calculation at point x. Sample is done by range search. 
         * {@see SampleOperation::sample_range}
         * 
         * @param feature space coordinate x (origin)
         * @param search parameters
         * @param kernel
         * @param weight function (defaults to NULL)
         * @param flag, indicating if the returned vector should be rounded
         *        to the coordinate system's resolution
         * @return mean shift vector
         */
        vector<T> 
        meanshift(const vector<T> &x,
                  const SearchParameters *params,
                  const Kernel<T> *kernel = new GaussianNormalKernel<T>(),
                  const WeightFunction<T> *w = NULL,
                  const bool normalize_shift = true);
    };
}

#endif
