#ifndef _M3D_Point_IMPL_H_
#define _M3D_Point_IMPL_H_

namespace m3D {

	using ::cfa::meanshift::Point;
    
    template <class T>
    M3DPoint<T>::M3DPoint() : Point<T>(), cluster(NULL)
    {
    }
    
    template <class T>
    M3DPoint<T>::M3DPoint( vector<T> &coord, vector<T>& value ) : Point<T>(coord,value), cluster(NULL)
    {
    }
    
    template <class T>
    M3DPoint<T>::M3DPoint( const M3DPoint<T> &o ) : Point<T>(o), cluster( o.cluster )
    {}
    
    template <class T>
    M3DPoint<T>
    M3DPoint<T>::operator = (const M3DPoint& o)
    {
    	M3DPoint<T> copy( o );
        
        return copy;
    }
    
    template <class T>
    M3DPoint<T>::~M3DPoint()
    {
    }
    
    template <class T>
    bool M3DPoint<T>::operator == (const M3DPoint<T> &o)
    {
        return this->values == o.values && cluster = o.cluster;
    }
    
};
    
#endif
