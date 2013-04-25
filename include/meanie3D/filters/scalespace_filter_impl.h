
#ifndef _M3D_ScaleSpaceFilter_Impl_H_
#define _M3D_ScaleSpaceFilter_Impl_H_

#include <exception>
#include <stdexcept>
#include <cmath>
#include <map>
#include <string>

#include <boost/progress.hpp>

#include <cf-algorithms/cf-algorithms.h>

#include <meanie3D/utils.h>

namespace m3D {

	using namespace std;
	using namespace ::cfa::meanshift;
	using ::cfa::utils::VisitUtils;
	using ::cfa::utils::coords::CoordinateSystem;

	using namespace std;
 
    template <typename T>
    ScaleSpaceFilter<T>::ScaleSpaceFilter(T scale,
                                          const vector<T> &resolution,
                                          T decay,
                                          bool show_progress)
    : FeatureSpaceFilter<T>(show_progress)
    , m_scale(scale), m_decay(decay), m_progress_bar(NULL)
    {
        if ( scale < 0 )
        {
            throw logic_error("scale can not be less than zero");
        }
        
        if ( decay < 0 || decay >= 1 )
        {
            throw logic_error("decay must be > 0 and < 1");
        }
        
        T filter_width = sqrt(ceil(-2.0*scale*log(decay)))/2;
        
        m_kernels.clear();
        
        for (size_t i = 0; i < resolution.size(); i++)
        {
            // calculate the distances vector
            
            size_t mask_size = filter_width / resolution[i];
            
            vector<T> distances(mask_size,0.0);
            
            for (size_t j=0; j < mask_size; j++)
            {
                distances[j] = ((T)j) * resolution[i];
            }
            
            // create the kernel
            
            ScaleSpaceKernel<T> kernel(scale,distances);
            
            m_kernels.push_back(kernel);
        }
    };

    template <typename T>
    ScaleSpaceFilter<T>::~ScaleSpaceFilter()
    {
    };
    
#pragma mark -
#pragma mark Abstract filter method
    
//    template <typename T>
//    void
//    ScaleSpaceFilter<T>::applyWithoutNewPoints( FeatureSpace<T> *fs )
//    {
//      using namespace cfa::utils::timer;
//      using namespace cfa::utils::vectors;
// 
//#if WRITE_FEATURESPACE
//        std::string fn = fs->filename() + "_scale_0.vtk";
//        VisitUtils<T>::write_featurespace_vtk( fn, fs );
//#endif
//        if ( m_scale > 0.0 )
//        {
//            // Replace the value of each point in feature-space with the average
//            // of the neighbours within a range determined by the scale parameter
//            
//            double filter_width = sqrt( ceil( -2.0 * m_scale * log(0.01) ) ) / 2;
//            
//            vector<T> spatial_range( fs->coordinate_system->size(), filter_width );
//            
//            boost::progress_display *progress_bar = NULL;
//            
//            if ( this->show_progress() )
//            {
//                cout << endl << "Applying scale filter t=" << m_scale << " (ranges " << spatial_range << ") ...";
//
//                progress_bar = new boost::progress_display( fs->size() );
//
//                start_timer();
//            }
//            
//            RangeSearchParams<T> params( spatial_range );
//            
//            // Make a copy of the feature space
//            
//            FeatureSpace<T> fs_copy = FeatureSpace<T>( *fs );
//            
//            // create a 'spatial only' index of the feature-space points
//            
//            PointIndex<T> *index = PointIndex<T>::create( &fs_copy, fs->coordinate_system->dimension_variables() );
//
//            // Initialize limit trackers
//            
//            map<int,T> min;
//            map<int,T> max;
//
//            // Pre-calculate the gauging factor for the gaussian kernel
//            
//            double gauging_factor = 1.0 / std::pow( (double)sqrt(2.0 * M_PI * m_scale), (double)fs->coordinate_system->size() );
//
//            // Iterate over feature-space
//
//            typename Point<T>::list::iterator pit;
//            
//            for ( pit = fs->points.begin(); pit != fs->points.end(); pit++ )
//            {
//                if ( this->show_progress() )
//                {
//                    progress_bar->operator++();
//                }
//                
//                Point<T> *p = *pit;
//                
//                // search in range filter_width around p in spatial index
//                
//                typename Point<T>::list *sample = index->search( p->coordinate, &params );
//                
//                // Iterate over the range variables in the sample
//                
//                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
//                {
//                    // replace the point's value at var_index with the gauss-weighted (by distance)
//                    // weighed sum of it's neighbours
//                    
//                    T sum = 0.0;
//                    
//                    typename Point<T>::list::iterator sample_iter;
//                    
//                    for ( sample_iter = sample->begin(); sample_iter != sample->end(); sample_iter++ )
//                    {
//                        typename Point<T>::ptr sample_point = *sample_iter;
//                        
//                        // distance from origin?
//                        
//                        T r = vector_norm( p->coordinate - sample_point->coordinate );
//                        
//                        double gauss = gauging_factor * exp(-(r*r)/(2*m_scale));
//                        
//                        // weighted sum
//                        
//                        sum += gauss * sample_point->values[var_index];
//                    }
//                    
//                    // replace the value by the kernel-weighted sum
//                    
//                    p->values[var_index] = sum;
//                    
//                    // make sure limit trackers are initialized
//                    
//                    typename map<int,T>::iterator m;
//                    
//                    if ( (m = min.find(var_index)) == min.end() )
//                    {
//                        min[var_index] = std::numeric_limits<T>::max();
//                    }
//                    
//                    if ( (m = max.find(var_index)) == max.end() )
//                    {
//                        max[var_index] = std::numeric_limits<T>::min();
//                    }
//                    
//                    // track limits
//                    
//                    if (p->values[var_index] < min[var_index] )
//                    {
//                        min[var_index] = p->values[var_index];
//                    }
//                    
//                    if (p->values[var_index] > max[var_index] )
//                    {
//                        max[var_index] = p->values[var_index];
//                    }
//                }
//                
//                delete sample;
//            }
//            
//            delete index;
//            
//            if ( this->show_progress() )
//            {
//                delete progress_bar;
//                cout << "done. (" << stop_timer() << "s)" << endl;
//                cout << "Value ranges of smoothed features:" << endl;
//                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
//                {
//                    cout << "\t#" << var_index << " : [" << min[var_index] << "," << max[var_index] << "]" << endl;
//                }
//            }
//        }
//        
//#if WRITE_FEATURESPACE
//        fn = fs->filename() + "_scale_" + boost::lexical_cast<string>(m_scale) + ".vtk";
//        VisitUtils<T>::write_featurespace_vtk( fn, fs );
//#endif
//    }
//    
//    
//    template <typename T>
//    void
//    ScaleSpaceFilter<T>::applyWithNewPointsRecursive( FeatureSpace<T> *fs,
//                                                      size_t dimensionIndex,
//                                                      typename cfa::utils::CoordinateSystem<T>::GridPoint &gridpoint )
//    {
//        using namespace std;
//        
//        NcDim dim = fs->coordinate_system->dimensions()[dimensionIndex];
//        
//        // iterate over dimensions
//        
//        for ( int index = 0; index < dim.getSize(); index++ )
//        {
//            // at each point for each variable
//
//            gridpoint[dimensionIndex] = index;
//            
//            if ( dimensionIndex < (gridpoint.size()-1) )
//            {
//                // recurse
//                
//                applyWithNewPointsRecursive(fs, dimensionIndex+1, gridpoint);
//            }
//            else
//            {
//                bool isOriginalPoint = false;
//                
//                if ( this->show_progress() )
//                {
//                    m_progress_bar->operator++();
//                }
//            
//                // Obtain the coordinate
//                
//                vector<T> x(gridpoint.size(),0);
//                
//                fs->coordinate_system->lookup(gridpoint,x);
//                
//                // search in range filter_width around p in spatial index
//                // Iterate over the range variables in the sample
//                // calculate the result from the points found in the index
//
//                typename Point<T>::list *sample = m_index->search( x, m_range );
//                
//                if ( sample->size() == 0 )
//                {
//                    delete sample;
//                    
//                    continue;
//                }
//                
//                // calculate smoothed value as gauss-weighed sum of
//                // the values in the sample around x
//
//                vector<T> values( fs->feature_variables().size(), 0.0 );
//                
//                typename Point<T>::list::iterator sample_iter;
//                
//                for ( sample_iter = sample->begin(); sample_iter != sample->end(); sample_iter++ )
//                {
//                    
//                    // calculate weight
//                    
//                    typename Point<T>::ptr sample_point = *sample_iter;
//                    
//                    // Remember, if a point existed at this coordinate
//                    // in the original feature-space
//                    
//                    if ( sample_point->coordinate == x )
//                    {
//                        isOriginalPoint = true;
//                    }
//
//                    vector<T> dx = x - sample_point->coordinate;
//                    
//                    T weight = m_gauging_factor * std::exp( - (dx*dx) / (2*m_scale) );
//                    
//                    // calculate values
//                    
//                    for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
//                    {
//                        values[var_index] += weight * sample_point->values[var_index];
//                    }
//                }
//                
//                // Track limits
//                
//                typename map<int,T>::iterator m;
//                
//                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
//                {
//                    if ( (m = m_min.find(var_index)) == m_min.end() )
//                    {
//                        m_min[var_index] = std::numeric_limits<T>::max();
//                    }
//                    
//                    if ( (m = m_max.find(var_index)) == m_max.end() )
//                    {
//                        m_max[var_index] = std::numeric_limits<T>::min();
//                    }
//                    
//                    if (values[var_index] < m_min[var_index] )
//                    {
//                        m_min[var_index] = values[var_index];
//                    }
//                    
//                    if (values[var_index] > m_max[var_index] )
//                    {
//                        m_max[var_index] = values[var_index];
//                    }
//                }
//                
//                delete sample;
//
//                // Create point
//                
//                for ( size_t i = 0; i < fs->coordinate_system->size(); i++ )
//                {
//                    values[i] = x[i];
//                }
//            
//                typename Point<T>::ptr p = PointFactory<T>::get_instance()->create(gridpoint,x,values);
//                
//                M3DPoint<T> *mp = (M3DPoint<T> *)p;
//                
//                mp->setIsOriginalPoint(isOriginalPoint);
//                
//                fs->points.push_back(p);
//            }
//        }
//    }
//    
//    
//    template <typename T>
//    void
//    ScaleSpaceFilter<T>::applyWithNewPoints( FeatureSpace<T> *fs )
//    {
//        using namespace cfa::utils::timer;
//        using namespace cfa::utils::vectors;
//        
//        using namespace cfa::utils::timer;
//        using namespace cfa::utils::vectors;
//        
//#if WRITE_FEATURESPACE
//        std::string fn = fs->filename() + "_scale_0.vtk";
//        VisitUtils<T>::write_featurespace_vtk( fn, fs );
//#endif
//        if ( m_scale > 0.0 )
//        {
//            // Replace the value of each point in feature-space with the average
//            // of the neighbours within a range determined by the scale parameter
//            
//            m_filter_width = sqrt( ceil( -2.0 * m_scale * log(0.01) ) ) / 2;
//            
//            // Pre-calculate the gauging factor for the gaussian kernel
//            
//            m_gauging_factor = 1.0 / std::pow( 2.0 * M_PI * m_scale, (double)fs->coordinate_system->size()/2 );
//
//            vector<T> spatial_range( fs->coordinate_system->size(), m_filter_width );
//            
//            m_range = new RangeSearchParams<T>( spatial_range );
//            
//            if ( this->show_progress() )
//            {
//                cout << endl << "Applying scale filter t=" << m_scale << " (ranges " << spatial_range << ") ...";
//                
//                long numPoints = 1;
//                
//                for ( size_t i=0; i < fs->coordinate_system->dimensions().size(); i++)
//                {
//                    numPoints *= fs->coordinate_system->dimensions()[i].getSize();
//                }
//                
//                m_progress_bar = new boost::progress_display( numPoints );
//                
//                start_timer();
//            }
//            
//            // create a 'spatial only' index of the feature-space points
//            
//            m_index = PointIndex<T>::create( fs, fs->coordinate_system->dimension_variables() );
//            
//            // Iterate over the dimensions
//            
//            typename cfa::utils::CoordinateSystem<T>::GridPoint gp = fs->coordinate_system->newGridPoint();
//            
//            FeatureSpace<T> *filtered_featurespace = new FeatureSpace<T>(*fs,false);
//            
//            this->applyWithNewPointsRecursive(filtered_featurespace, 0, gp);
//            
//            // replace points in the featurespace with
//            // the filtered results
//            
//            fs->clear_points();
//            
//            fs->points = filtered_featurespace->points;
//            
//            // Finish 
//
//            if ( this->show_progress() )
//            {
//                cout << "done. (" << stop_timer() << "s)" << endl;
//                cout << "Value ranges of smoothed features:" << endl;
//                for ( size_t var_index = fs->coordinate_system->size(); var_index < fs->feature_variables().size(); var_index++ )
//                {
//                    cout << "\t#" << var_index << " : [" << m_min[var_index] << "," << m_max[var_index] << "]" << endl;
//                }
//
//                delete m_progress_bar;
//                m_progress_bar = NULL;
//            }
//#if WRITE_FEATURESPACE
//            fn = fs->filename() + "_scale_" + boost::lexical_cast<string>(m_scale) + ".vtk";
//            VisitUtils<T>::write_featurespace_vtk( fn, fs );
//#endif
//        }
//    }
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndexRecursive(FeatureSpace<T> *fs,
                                                      ArrayIndex<T> *originalIndex,
                                                      ArrayIndex<T> *filteredPoints,
                                                      vector<size_t> dimensionIndexes,
                                                      size_t fixedDimensionIndex,
                                                      size_t dimensionIndex,
                                                      typename CoordinateSystem<T>::GridPoint& gridpoint)
    {
        // loop over the dimension with the given index. For each
        // coordinate, loop over the other two. At Each point, apply
        // the one-dimensional kernel

        using namespace std;
        
        CoordinateSystem<T> *cs = fs->coordinate_system;
        
        size_t realDimIndex = dimensionIndexes[dimensionIndex];
        
        NcDim dim = cs->dimensions()[realDimIndex];

        // iterate over dimensions

        for ( int index = 0; index < dim.getSize(); index++ )
        {
            // at each point for each variable

            gridpoint[realDimIndex] = index;

            if ( dimensionIndex < (gridpoint.size()-1) )
            {
                // recurse

                applyWithArrayIndexRecursive(fs, originalIndex, filteredPoints, dimensionIndexes, fixedDimensionIndex, dimensionIndex+1, gridpoint);
            }
            else
            {
                // we reached the fixed dimension!
                
                if ( this->show_progress() )
                {
                    m_progress_bar->operator++();
                }
                
                ScaleSpaceKernel<T> g = this->m_kernels[realDimIndex];
                
                // Find the boundaries. Take care not to step
                // outside the bounds of the array
                
                int width = g.values().size() - 1;
                
                int gpIndex = (int)gridpoint[realDimIndex];
                
                int minIndex = (gpIndex - width >= 0) ? (gpIndex - width) : 0;
                
                int maxIndex = ((gpIndex + width) < (dim.getSize()-1)) ? (gpIndex + width) : (dim.getSize()-1);
                
                // Convolute in 1D around the given point with
                // the mask size determined by the kernel
                // Run the convolution for each feature variable
                
                typename CoordinateSystem<T>::GridPoint gridIter = gridpoint;
                
                for (int varIndex=0; varIndex < fs->variables().size(); varIndex++)
                {
                    T sum = 0.0;
                    
                    size_t sumCount = 0;
                    
                    for (int i=minIndex; i<maxIndex; i++)
                    {
                        // set gridpoint to current position
                        
                        gridIter[realDimIndex] = i;
                        
                        // get the point at the current position from
                        // the array index

                        Point<T> *pIter = originalIndex->get(gridIter);
                        
                        if ( pIter == NULL ) continue;
                        
                        // apply the pre-sampled gaussian and sum up

                        size_t d = (i <= index) ? (index-i) : (i-index);
                        
                        sum += g.value(d) * pIter->values[cs->size()+varIndex];
                        
                        sumCount++;
                    }
                    
                    // No muss, no fuss?
                    
                    if (sumCount==0) continue;
                    
                    // Fuss! Fetch the point from the array index
                    
                    Point<T> *p = filteredPoints->get(gridpoint);
                    
                    // If no point existed, decide if we need to create one
                    
                    if ( p == NULL )
                    {
                        // Create a new point with default values
                        // and insert into array index
                        
                        typename CoordinateSystem<T>::Coordinate coordinate = cs->newCoordinate();
                        
                        cs->lookup(gridpoint,coordinate);
                        
                        vector<T> values = coordinate;
                        
                        values.resize(fs->feature_variables().size(),0.0);
                        
                        p = PointFactory<T>::get_instance()->create(gridpoint,coordinate,values);
                        
                        filteredPoints->set(gridpoint,p);
                        
                        m_created_points++;
                    }
                    else
                    {
                        m_modified_points++;
                    }
                    
                    // If we have a point after all that, update it with the
                    // filtered value
                    
                    if ( p != NULL )
                    {
                        p->values[fs->coordinate_system->size()+varIndex] = sum;
                    }
                }
            }
        }
    }
            
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndexForDimension(FeatureSpace<T> *fs,
                                                         ArrayIndex<T> *originalIndex,
                                                         ArrayIndex<T> *filteredPoints,
                                                         size_t fixedDimensionIndex)
    {
        // iterate over the other two
        
        // Create a vector, that enumerates the dimensions such, that the
        // fixed dimension comes last
     
        vector<size_t> dimensionIndexes;
        
        for (size_t j=0; j<fs->coordinate_system->size(); j++)
        {
            if (j==fixedDimensionIndex) continue;
            dimensionIndexes.push_back(j);
        }
        
        dimensionIndexes.push_back(fixedDimensionIndex);
        
        typename CoordinateSystem<T>::GridPoint gridpoint = fs->coordinate_system->newGridPoint();
        
        // Now recurse into the structure, bearing in mind that
        // the dimensions have been re-ordered
        
        applyWithArrayIndexRecursive(fs, originalIndex, filteredPoints, dimensionIndexes, fixedDimensionIndex, 0, gridpoint);
    }
    
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndex(FeatureSpace<T> *fs)
    {
        using namespace std;
        
        // index the original
        
        if ( this->show_progress() )
        {
            cout << "Constructing array index ...";
            
            start_timer();
        }

        ArrayIndex<T> *originalIndex = new ArrayIndex<T>(fs);
        
        if ( this->show_progress() )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;
        }
        
        // Create a copy and index that as well. The results
        // are going to be stored in the index
        
        if ( this->show_progress() )
        {
            cout << "Copying feature-space ...";
            
            start_timer();
        }
        
        FeatureSpace<T> *filteredSpace = new FeatureSpace<T>(fs,true);
        
        if ( this->show_progress() )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;
        }
        
        if ( this->show_progress() )
        {
            cout << "Constructing array index ...";
            
            start_timer();
        }

        ArrayIndex<T> *filteredIndex = new ArrayIndex<T>(filteredSpace);
        
        if ( this->show_progress() )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;
        }

        if ( this->show_progress() )
        {
            cout << endl << "Applying scale filter t=" << m_scale << " decay=" << m_decay << " ... " << endl;

            long numPoints = 1;

            for ( size_t i=0; i < fs->coordinate_system->dimensions().size(); i++)
            {
                numPoints *= fs->coordinate_system->dimensions()[i].getSize();
            }

            m_progress_bar = new boost::progress_display( fs->coordinate_system->dimensions().size() * numPoints );
            
            start_timer();
        }
        
        m_modified_points = m_created_points = 0;
        
        for (size_t dimIndex=0; dimIndex < fs->coordinate_system->size(); dimIndex++)
        {
            applyWithArrayIndexForDimension(fs, originalIndex, filteredIndex, dimIndex);
            
            if (dimIndex < fs->coordinate_system->size() - 1)
            {
                delete originalIndex;
            
                originalIndex = new ArrayIndex<T>(filteredIndex);
            }
        }
        
        // replace the points in the original with the filtered
        // array index results
        
        filteredIndex->replace_points(fs);
        
        if ( this->show_progress() )
        {
            cout << "done. (" << stop_timer() << "s)" << endl;
            cout << "Filtered featurespace contains " << fs->size() << " points."
                 << "(" << m_created_points << " new points)" << endl;
            delete m_progress_bar;
            m_progress_bar = NULL;
        }
#if WRITE_FEATURESPACE
        string fn = fs->filename() + "_scale_" + boost::lexical_cast<string>(m_scale) + ".vtk";
        VisitUtils<T>::write_featurespace_vtk( fn, fs );
#endif
        
        // Clean up
        
        delete filteredIndex;
        
        delete filteredSpace;
    }
    
    template <typename T>
    void
    ScaleSpaceFilter<T>::apply( FeatureSpace<T> *fs )
    {
        this->applyWithArrayIndex(fs);
    }
    
};

#endif
