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
        
        virtual Point<T> * create() {
            return new M3DPoint<T>();
        }
        
        virtual Point<T> * create( vector<size_t> &gridpoint, vector<T> &coord, vector<T> &value )
        {
            M3DPoint<T> *p = new M3DPoint<T>(gridpoint,coord,value);
            
            return p;
        }
        
        virtual Point<T> * create( vector<T> &coord, vector<T> &value )
        {
            M3DPoint<T> *p = new M3DPoint<T>(coord,value);
            
            return p;
        }

        virtual Point<T> *
        copy( const Point<T> *p)
        {
            return new M3DPoint<T>( (M3DPoint<T>*) p );
        }

    };
    
}

#endif
