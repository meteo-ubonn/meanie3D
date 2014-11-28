/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


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
    
    template <class T>
    void 
    Point<T>::print(unsigned short num_tabs)
    {
        for (unsigned short ti=0; ti<num_tabs; ti++)
            cout << "\t";
        
        std::cout << "gridpoint=" << this->gridpoint
                << " coordinate=" << this->coordinate
                << " values=" << this->values
                << endl;
    }

}

#endif
