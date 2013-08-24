#ifndef _M3D_Point_H_
#define _M3D_Point_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <cf-algorithms/featurespace/point.h>

namespace m3D {

	using cfa::meanshift::Point;

	// Forward Declaration

	template <typename T>
	class Cluster;

	/** This represents one point f in feature space F.
     */
    template <class T> 
    struct M3DPoint : public Point<T>
    {
        
    private:
        

    public:
        
#pragma mark -
#pragma mark Definitions
        
        typedef M3DPoint<T> * ptr;
        typedef vector< ptr > list;

#pragma mark -
#pragma mark public properties
        
        typename Cluster<T>::ptr    cluster;
        bool                        isBoundary;
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Default constructor
         */
    	M3DPoint()
        : Point<T>()
        , cluster(NULL)
        , isBoundary(false)
        {};

        /** Constructor.
         * @param gridpoint
         * @param coordinate
         * @param values
         */
        M3DPoint( vector<size_t> gridpoint, vector<T> &coord, vector<T>& value )
        : Point<T>(gridpoint,coord,value)
        , cluster(NULL)
        , isBoundary(false)
        {};

        /** Constructor.
         * @param gridpoint
         * @param coordinate
         * @deprecated
         */
        M3DPoint( vector<T> &coord, vector<T>& value )
        : Point<T>(coord,value)
        , cluster(NULL)
        , isBoundary(false)
        {};

        /** Copy constructor
         */
        M3DPoint( const M3DPoint<T> &o )
        : Point<T>(o)
        , cluster( o.cluster )
        , isBoundary(o.isBoundary)
        {};
        
        /** Copy constructor on pointer
         */
        M3DPoint( const M3DPoint<T> *o )
        : Point<T>(o)
        , cluster( o->cluster )
        , isBoundary(o->isBoundary)
        {};

        /** Copy operator */
        
        M3DPoint<T>
        operator = (const M3DPoint& o) { M3DPoint<T> copy( o ); return copy; }

        /** Destructor */
        virtual ~M3DPoint()
        {
            //cout << "~MDPoint(" << this << ")" << endl;
        }
        
#pragma mark -
#pragma mark Operators
        
        /** The equal operator is overloaded. A data point is equal to another data
         * point if it has the same coordinate. This is necessary to be able to get 
         * rid of the crappy stl <map> implementation with it's gazillion calls to
         * copy constructors etc.
         */
        bool
        operator == (const M3DPoint<T> &o) { return this->values == o.values; };
        
    };
}

#endif
