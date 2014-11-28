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

#ifndef M3D_HISTOGRAM_H
#define M3D_HISTOGRAM_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/featurespace.h>

#include <stdlib.h>
#include <vector>

namespace m3D {

    /** This represents one point f in feature space F.
     */
    template <class T> 
    struct Histogram
    {
    private:

        vector<size_t>  m_bins;

        bool            m_sum_dirty;

        size_t          m_sum;

    public:

        typedef Histogram<T>* ptr;

#pragma mark -
#pragma mark Constructor/Destructor

        /** Default @constructor is private
         */
        Histogram() : m_sum_dirty(true),m_sum(0) {};

        /** @constructor
         * @param number of bins
         */
        Histogram( const size_t &size ) : m_bins(vector<size_t>(size,0)), m_sum_dirty(true), m_sum(0)  {};

        /** Constructor.
         * @param initial bins
         */
        Histogram( vector<size_t> &bins ) : m_bins(bins), m_sum_dirty(true), m_sum(0) {};

        /** Copy constructor
         */
        Histogram( const Histogram<T> &o ) : m_bins(o.bins()), m_sum_dirty(true), m_sum(0) {};

        /** Destructor 
         */
        ~Histogram()
        {
            cout << "~Histogram()" << endl;
        };

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
        const size_t size() const { return this->m_bins.size(); };

        /** Contains the sum of all bins 
         */
        const size_t sum();

        /** Access all bins at once
         */
        const vector<size_t> &bins() const { return m_bins; };

        /** Access a bin directly
         * @param index
         */
        T& operator [] (const size_t index);

        /** Allocates an int array and fills it with the bins. 
         * @return int array of size this->size()
         */
        int* bins_as_int_array() const;

        /** Allocates a float array and fills it with the bins.
         * @return float array of size this->size()
         */
        float* bins_as_float_array() const;

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
        static typename Histogram<T>::ptr
        create( typename Point<T>::list &points, size_t variable_index, T min, T max, size_t number_of_bins = 10 );

#pragma mark -
#pragma mark Histogram Correlation

        /** Spearman histogram corellation with another histogram
         * (see numerical recipes 2nd edition)
         * @param other histogram
         * @return rho [-1.0 .. 1.0]
         */
        T correlate_spearman( const typename Histogram<T>::ptr o );

        /** Kendall histogram corellation with another histogram
         * (see numerical recipes 2nd edition)
         * @param other histogram
         * @return tau [-1.0 .. 1.0]
         */
        T correlate_kendall( const typename Histogram<T>::ptr o );

    };
}

#endif
