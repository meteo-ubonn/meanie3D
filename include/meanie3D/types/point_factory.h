#ifndef _M3D_PointFactory_H_
#define _M3D_PointFactory_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <meanie3D/types/point.h>
#include <cf-algorithms/types.h>
#include <cf-algorithms/featurespace/point_factory.h>

namespace m3D {
    
	using  ::cfa::meanshift::PointFactory;

    /** Default point factory
     */
    template <typename T>
    class M3DPointFactory : public PointFactory<T>
    {
        
        /** This is an abstract method used to create new instances of Point.
         */
        Point<T> * create() {
            return new M3DPoint<T>();
        }
        
        /** This is an abstract method used to create new instances of Point.
         * @param vector with initial point coordinates
         * @param vector with initial point value
         */
        Point<T> * create( vector<T> &coord, vector<T> &value ) {
            return new M3DPoint<T>(coord,value);
        }
    };
    
};

#endif
