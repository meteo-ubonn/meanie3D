#ifndef M3D_WEIGHTFUNCTION_H
#define M3D_WEIGHTFUNCTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>

namespace m3D { 

    /** Weight function interface. The weight function plays an important
     * role in the mean-shift algorithm. It replaces density as criterion
     * in places, where the actual spatial density is homogenous, as is
     * often the case in scientific, gridded data sets.
     */
    template <typename T>
    class WeightFunction {
    public:

        /** Weight at given point in feature-space 
         * @param point object
         * @return weight
         */
        virtual T operator()(const typename Point<T>::ptr p) const = 0;

        virtual ~WeightFunction() {}
    };
}

#endif
