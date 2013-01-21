#ifndef _M3D_Histogram_H_
#define _M3D_Histogram_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <stdlib.h>
#include <vector>

#include <cf-algorithms/featurespace/point.h>

namespace m3D {
    
    using cfa::meanshift::Point;

	/** This represents one point f in feature space F.
     */
    template <class T> 
    struct Histogram
    {
    public:

        /** Different ways of comparing histograms 
         */
        typedef enum {
            HistogramCorrelationSpearman
        } HistogramCorrelation;
        
    private:
        
        vector<size_t>   m_bins;
        
        /** Default @constructor is private
         */
    	Histogram() {};
        
    public:
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** @constructor
         * @param number of bins
         */
    	Histogram( const size_t &size ) : m_bins(vector<size_t>(size,0)) {};

        /** Constructor.
         * @param initial bins
         */
        Histogram( vector<size_t> &bins ) {this->m_bins = bins; };

        /** Copy constructor
         */
        Histogram( Histogram<T> &o ) : m_bins(o.bins()) {};
        
        /** Destructor 
         */
        ~Histogram() {};
        
#pragma mark -
#pragma mark Operators

        /** equals operator == */
        
        bool
        operator == (const Histogram<T> &o) { return this->m_bins = o.bins(); };
        
        /** assignment operator = */
        
        Histogram<T>
        operator = (const Histogram& o) { return Histogram<T>(o); };

#pragma mark -
#pragma mark Accessors
        
        /** Get the number of bins
         * @return size
         */
        const size_t size() { return this->m_bins.size(); };
        
        /** Access all bins at once
         */
        vector<size_t> &bins() { return m_bins; };
        
        /** Access a bin directly
         * @param index
         */
        T& operator [] (const size_t index) { return this->m_bins[index]; }
        
#pragma mark -
#pragma mark Comparisons
        
        /** Correlates this histogram with another. 
         * @param histogram
         * @param correlation method
         * @return correlation factor [0..1].
         */
        T correlate( const Histogram<T> &other, HistogramCorrelation method = HistogramCorrelationSpearman );
        
#pragma mark -
#pragma mark Factory Methods
        
        /** Creates a histogram from a point list. Classes are created
         * equidistanty in intervals of size (max-min)/number_of_classes
         * @param list
         * @param which variable should be indexed
         * @param lowest value in the histogram classes
         * @param highest value in the histogram classes
         * @param number of bins (default 10).
         */
        static Histogram<T> *
        create( typename Point<T>::list &points, size_t variable_index, T min, T max, size_t number_of_bins = 10 );
        
    };
};

#endif
