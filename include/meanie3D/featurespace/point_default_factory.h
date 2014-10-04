#ifndef M3D_POINT_DEFAULT_FACTORY_H
#define M3D_POINT_DEFAULT_FACTORY_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D {

    // Forward declaration
    template <typename T>
    class PointFactory;

    /** Default point factory
     */
    template <typename T>
    class PointDefaultFactory : public PointFactory<T>
    {

        virtual Point<T> * create() {
            return new Point<T>();
        }

        virtual Point<T> * create( vector<int> &gridpoint, vector<T> &coord, vector<T> &value ) {
            return new Point<T>(gridpoint,coord,value);
        }

        virtual Point<T> * create( vector<T> &coord, vector<T> &value ) {
            return new Point<T>(coord,value);
        }

        virtual Point<T> *
        copy( const Point<T> *p)
        {
            return new Point<T>(p);
        }

    };
}

#endif
