#ifndef _M3D_Histogram_Impl_H_
#define _M3D_Histogram_Impl_H_

#include <vector>
#include <utility>

#include <numericalrecipes/nr.h>

#include <cf-algorithms/featurespace/point.h>

extern "C"
{
    extern void spear(float data1[],float data2[],unsigned long n, float*d,float*zd, float*probd,float*rs,float*probrs);
    extern void kendl1(float data1[], float data2[], unsigned long n, float *tau, float *z, float *prob);
}


namespace m3D {

	using ::cfa::meanshift::Point;
    using namespace std;
    
    template <typename T>
    int *
    Histogram<T>::bins_as_int_array() const
    {
        int *array = new int( this->size() );
        
        for (size_t i=0; i < this->m_bins.size(); i++ )
        {
            array[i] = this->m_bins[i];
        }
        
        return array;
    };
    
    template <typename T>
    float *
    Histogram<T>::bins_as_float_array() const
    {
        float *array = (float *) malloc( sizeof(float) * this->m_bins.size() );
        
        for (size_t i=0; i < this->m_bins.size(); i++ )
        {
            array[i] = (float)this->m_bins[i];
        }
        
        return array;
    };

    template <typename T>
    T
    Histogram<T>::correlate_spearman( const typename Histogram<T>::ptr o )
    {
        // Transfer data into correct typed arrays
        
        float *h1 = this->bins_as_float_array();
        
        float *h2 = o->bins_as_float_array();
        
        float d,zd,probd,probrs,rho;
        
        spear( h1, h2, this->size(), &d, &zd, &probd, &rho, &probrs );
        
#if DEBUG_TRACKING
        cout << "spearman rank correlation of " << this->bins() << " with " << o->bins();
        cout << "  =>  rho=" << rho << " d=" << d << " zd=" << zd << " probrs=" << probrs << endl;
#endif
        return isnan(rho) ? (T)-1.0 : (T)rho;
    }
    
    template <typename T>
    T
    Histogram<T>::correlate_kendall( const typename Histogram<T>::ptr o )
    {
        // Transfer data into correct typed arrays
        
        float *h1 = this->bins_as_float_array();
        
        float *h2 = o->bins_as_float_array();
        
        float z,prob,tau;
        
        kendl1( h1, h2, (int)this->size(), &tau, &z, &prob );

#if DEBUG_TRACKING
        cout << "kendall  rank correlation of " << this->bins() << " with " << o->bins();
        cout << "  =>  tau=" << tau << " z=" << z << " prob=" << prob << endl;
#endif
        return isnan(tau) ? (T)-1.0 : (T)tau;
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
    typename Histogram<T>::ptr
    Histogram<T>::create( typename Point<T>::list &points, size_t variable_index, T min, T max, size_t number_of_bins )
    {
        if ( max == min )
        {
            cerr << "Histogram::create:ERROR:degenerate case, min==max" << endl;
            
            vector<size_t> bins(1,points.size());
            
            return new Histogram<T>(bins);
        }
        
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
        
        return new Histogram<T>( bins );
    
    };

};

#endif
