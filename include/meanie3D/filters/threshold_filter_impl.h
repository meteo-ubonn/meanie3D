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
