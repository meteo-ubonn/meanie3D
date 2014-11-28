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
#ifndef M3D_THRESHOLDFILTER_IMPL_H
#define M3D_THRESHOLDFILTER_IMPL_H

#include <meanie3D/utils.h>

#include "threshold_filter.h"

namespace m3D {

    template <typename T>
    ThresholdFilter<T>::ThresholdFilter(const vector<T> &thresholds) : m_thresholds(thresholds) {}

    template <typename T>
    ThresholdFilter<T>::~ThresholdFilter() {}

#pragma mark -
#pragma mark Accessors

    template <typename T>
    void
    ThresholdFilter<T>::set_thresholds( const vector<T> &t ) { m_thresholds = t; }

    template <typename T>
    vector<T>
    ThresholdFilter<T>::thresholds() { return m_thresholds; }

#pragma mark -
#pragma mark Abstract filter method

    template <typename T>
    void
    ThresholdFilter<T>::apply( FeatureSpace<T> *fs )
    {
        using namespace std;

        // Check the dimensions
        assert( fs->dimension == m_thresholds.size() );

        boost::progress_display *progress_bar = NULL;

        if ( this->show_progress() )
        {
            cout << endl << "Applying threshold filter t=" << m_thresholds << " ...";

            progress_bar = new boost::progress_display( fs->size() );

            start_timer();
        }

        typename Point<T>::list::iterator pit;

        for ( pit = fs->points.begin(); pit != fs->points.end(); )
        {
            if ( this->show_progress() )
            {
                progress_bar->operator++();
            }

            Point<T> *p = *pit;

            bool threshold_ok = true;

            for ( size_t i=0; i < fs->dimemsion && threshold_ok; i++)
            {
                threshold_ok = p->values[i] >= m_thresholds[i];
            }

            if ( !threshold_ok )
            {
                fs->points.erase( pit );
            }
            else
            {
                pit++;
            }
        }

        if ( this->show_progress() )
        {
            delete progress_bar;
            cout << "done. (" << stop_timer() << "s)" << std::endl;
        }
    }
}

#endif
