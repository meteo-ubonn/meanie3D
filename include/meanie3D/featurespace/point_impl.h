#ifndef M3D_FEATURESPACEPOINT_IMPL_H
#define M3D_FEATURESPACEPOINT_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <iostream>
#include <math.h>

#include "point.h"

namespace m3D { 

    template <class T>
    Point<T>::Point()
    : trajectory_length(0)
    , isOriginalPoint(false)
    {
    }

    template <class T>
    Point<T>::Point( vector<int> &gp, vector<T> &coord, vector<T>& value )
    : coordinate(coord)
    , gridpoint(gp)
    , values(value)
    , isOriginalPoint(false)
    , cluster(NULL)
    , isBoundary(false)
    {
        trajectory_length = 0;
    }

    template <class T>
    Point<T>::Point( vector<T> &coord, vector<T>& value )
    : coordinate(coord)
    , values(value)
    , isOriginalPoint(false)
    , cluster(NULL)
    , isBoundary(false)
    {
        trajectory_length = 0;
    }

    template <class T>
    Point<T>::Point( const Point<T> &o )
    : coordinate( o.coordinate )
    , values( o.values )
    , trajectory_length( o.trajectory_length )
    , shift( o.shift )
    , gridded_shift(o.gridded_shift)
    , gridpoint(o.gridpoint)
    , isOriginalPoint(o.isOriginalPoint)
    , cluster(o.cluster)
    , isBoundary(o.isBoundary)
    {}

    template <class T>
    Point<T>::Point( const Point<T> *o )
    : coordinate( o->coordinate )
    , gridpoint(o->gridpoint)
    , values( o->values )
    , trajectory_length( o->trajectory_length )
    , shift( o->shift )
    , gridded_shift(o->gridded_shift)
    , isOriginalPoint(o->isOriginalPoint)
    , cluster(o->cluster)
    , isBoundary(o->isBoundary)
    {}

    template <class T>
    Point<T>
    Point<T>::operator = (const Point& o)
    {
        Point<T> copy( o );

        return copy;
    }

    template <class T>
    bool Point<T>::operator == (const Point<T> &o)
    {
        return values == o.values;
    }
}

#endif
