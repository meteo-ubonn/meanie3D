#ifndef _M3D_Histogram_Impl_H_
#define _M3D_Histogram_Impl_H_

#include <vector>
#include <utility>

#include <cf-algorithms/featurespace/point.h>

namespace m3D {

	using ::cfa::meanshift::Point;
    using namespace std;
    
    template <typename T>
    T correlate_spearman(const vector<T> v1, const vector<T> v2)
    {
        T result = 0.0;
        
        return result;
    }
    
    
    template <typename T>
    T
    Histogram<T>::correlate( const Histogram<T> &other, HistogramCorrelation method )
    {
        T result = 0;
        
        switch (method)
        {
            case HistogramCorrelationSpearman:
            {
                result = correlate_spearman( this->values(), other.values() );
            }
            break;
        }
        
        return result;
    }
    
#pragma mark -
#pragma mark Factory Methods
    
    /** Creates a histogram from a point list. Classes are created
     * equidistanty in intervals of size (max-min)/number_of_classes
     * @param list
     * @param which variable should be indexed
     * @param lowest value in the histogram classes
     * @param highest value in the histogram classes
     * @param number of classes.
     */
    template <typename T>
    Histogram<T> *
    Histogram<T>::create( typename Point<T>::list &points, size_t variable_index, T min, T max, size_t number_of_bins )
    {
        assert( max > min );
        
        typedef pair<T,T> class_t;
        
        typedef vector< class_t > classes_t;
        
        // create classes
        
        classes_t classes(number_of_bins);
        
        T dx = (max - min) / ((T) number_of_bins);
        
        for (size_t n = 0; n < number_of_bins; n++)
        {
            classes[n].first = min + n * dx;

            classes[n].second = min + (n+1) * dx;
        }
        
        // Count the number of points in each class
        
        vector<size_t> bins(number_of_bins,0);
        
        typename Point<T>::list::iterator it;
        
        for ( it = points.begin(); it != points.end(); it++ )
        {
            typename Point<T>::ptr p = *it;
            
            T value = p->values[variable_index];
            
            // TODO: obtain the class index through calculation
            
            if ( value >= max )
            {
                bins[ number_of_bins-1 ] += 1;
            }
            else
            {
                for ( size_t class_index = 0; class_index < classes.size(); class_index++ )
                {
                    class_t c = classes[class_index];
                    
                    if ( value >= c.first && value < c.second )
                    {
                        bins[class_index] += 1;
                        
                        break;
                    }
                }
            }
        }
        
        Histogram<T> *result = new Histogram<T>( bins );
    
        return result;
    
    };

};

#endif
