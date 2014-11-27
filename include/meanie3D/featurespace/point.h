#ifndef M3D_POINT_H
#define M3D_POINT_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering.h>
#include <meanie3D/parallel.h>
#include <meanie3D/utils.h>

#include <vector>

namespace m3D {

    using namespace std;
    
    // Forward Declaration
    template <typename T> class Cluster;

    /** This represents one point f in feature space F.
     */
    template <class T> 
    struct Point
    {

    private:

    public:

#pragma mark -
#pragma mark Definitions

        typedef Point<T> * ptr;
        typedef vector< ptr > list;
        typedef vector< ptr > set;

#pragma mark -
#pragma mark public properties

        /** spatial coordinate */
        vector<T>       coordinate;

        /** added for quicker array indexing */
        vector<int>     gridpoint;

        /** actual point in feature space */
        vector<T>       values;

        /** @deprecated */
        size_t          trajectory_length;

        /** meanshift vector (rounded to resolution) */
        vector<T>       shift;

        /** meanshift vector in grid points */
        vector<int>     gridded_shift;

        /** Remembers if this point was part of the original
         * feature-space construction
         */
        bool            isOriginalPoint;

        /** A pointer to the cluster this point belongs to.
         */
        Cluster<T>      *cluster;

        /** Flag indicating if this point is on the cluster's
         * boundary or not.
         */
        bool            isBoundary;

#if PROVIDE_MUTEX
        boost::mutex mutex;
#endif

#pragma mark -
#pragma mark Constructor/Destructor

        /** Copy constructor
         */
        Point( const Point<T> &o );

        /** Copy constructor on pointer
         */
        Point( const Point<T> *o );

        /** Copy operator */

        Point<T>
        operator = (const Point& o);

        /** Default constructor
         */
        Point();

        /** Constructor.
         * @param gridpoint
         * @param coordinate
         * @param values
         */
        Point(vector<int> &gridpoint, 
               vector<T> &coordinate, 
               vector<T>& values);

        /** Constructor.
         * @param coordinate
         * @param values
         * @deprecated
         */
        Point(vector<T> &coordinate,
              vector<T>& values );

        /** Destructor 
         */
        virtual ~Point() {}

#pragma mark -
#pragma mark Operators

        /** The equal operator is overloaded. A data point is equal to another data
         * point if it has the same coordinate. This is necessary to be able to get 
         * rid of the crappy stl <map> implementation with it's gazillion calls to
         * copy constructors etc.
         */
        bool
        operator == (const Point<T> &o);
        
#pragma mark -
#pragma mark Misc
        
        /** prints point info out to cout
         * @param number of tabs to prepend
         */
        void print(unsigned short num_tabs=0);
        

    };
}

#endif
