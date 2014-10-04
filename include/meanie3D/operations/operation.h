#ifndef M3D_OPERATION_H
#define M3D_OPERATION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace.h>
#include <meanie3D/index.h>

namespace m3D { 

    /** Abstract base class for Operations on FeatureSpace
     */
    template <typename T> 
    class Operation
    {
    protected:

        FeatureSpace<T>  *feature_space;

        PointIndex<T>    *point_index;

    public:

        Operation( FeatureSpace<T> *fs, PointIndex<T> *index ) : feature_space(fs),point_index(index) {}

        virtual ~Operation() {};
    };   
}
    
#endif
