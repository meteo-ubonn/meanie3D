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


#ifndef M3D_HISTOGRAM_IMPL_H
#define M3D_HISTOGRAM_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/numericalrecipes.h>

#include <vector>
#include <utility>

#include "histogram.h"

namespace m3D {

    using namespace std;

    template <typename T>
    const size_t
    Histogram<T>::sum()
    {
        size_t sum = 0;
        for (size_t i = 0; i < m_bins.size(); i++) {
            sum += m_bins[i];
        }
        return sum;
    }

    template <typename T>
    T&
    Histogram<T>::operator[](const size_t index)
    {
        return this->m_bins[index];
    }

    template <typename T>
    int *
    Histogram<T>::bins_as_int_array() const
    {
        int *array = new int( this->size());
        for (size_t i = 0; i < this->m_bins.size(); i++) {
            array[i] = this->m_bins[i];
        }
        return array;
    }

    template <typename T>
    float *
    Histogram<T>::bins_as_float_array() const
    {
        float *array = (float *) malloc(sizeof (float) * this->m_bins.size());
        for (size_t i = 0; i < this->m_bins.size(); i++) {
            array[i] = (float) this->m_bins[i];
        }
        return array;
    }

    template <typename T>
    T
    Histogram<T>::correlate_spearman(const typename Histogram<T>::ptr o)
    {
        // Transfer data into correct typed arrays
        float *h1 = this->bins_as_float_array();
        float *h2 = o->bins_as_float_array();
        float d, zd, probd, probrs, rho;
        spear(h1, h2, this->size(), &d, &zd, &probd, &rho, &probrs);
        return isnan(rho) ? (T) - 1.0 : (T) rho;
    }

    template <typename T>
    T
    Histogram<T>::correlate_kendall(const typename Histogram<T>::ptr o)
    {
        // Transfer data into correct typed arrays
        float *h1 = this->bins_as_float_array();
        float *h2 = o->bins_as_float_array();
        float z, prob, tau;
        kendl1(h1, h2, (int) this->size(), &tau, &z, &prob);
        return isnan(tau) ? (T) - 1.0 : (T) tau;
    }

#pragma mark -
#pragma mark Factory Methods

    /** 
     * Creates a histogram from a point list. Classes are created
     * equidistantly in intervals of size (max-min)/number_of_classes
     * 
     * @param list
     * @param which variable should be indexed
     * @param lowest value in the histogram classes
     * @param highest value in the histogram classes
     * @param number of classes.
     */
    template <typename T>
    typename Histogram<T>::ptr
    Histogram<T>::create(typename Point<T>::list &points, size_t variable_index, T min, T max, size_t number_of_bins)
    {
        if (max == min) {
            cerr << "ERROR:histogram::create:ERROR:degenerate case, min==max" << endl;
            vector<size_t> bins(1, points.size());
            return new Histogram<T>(bins);
        }

        typedef pair<T, T> class_t;
        typedef vector< class_t > classes_t;

        // create classes
        classes_t classes(number_of_bins);
        T dx = (max - min) / ((T) number_of_bins);
        for (size_t n = 0; n < number_of_bins; n++) {
            classes[n].first = min + n * dx;
            classes[n].second = min + (n + 1) * dx;
        }

        // Count the number of points in each class
        vector<size_t> bins(number_of_bins, 0);
        typename Point<T>::list::iterator it;
        for (it = points.begin(); it != points.end(); it++) {
            typename Point<T>::ptr p = *it;
            T value = p->values[variable_index];
            // TODO: obtain the class index through calculation
            if (value >= max) {
                bins[ number_of_bins - 1 ] += 1;
            } else {
                for (size_t class_index = 0;
                        class_index < classes.size(); class_index++) {
                    class_t c = classes[class_index];
                    if (value >= c.first && value < c.second) {
                        bins[class_index] = bins[class_index] + 1;
                        break;
                    }
                }
            }
        }
        return new Histogram<T>(bins);
    }
}

#endif
