#ifndef M3D_OPERATION_MEANSHIFTOPERATION_IMPL_H
#define M3D_OPERATION_MEANSHIFTOPERATION_IMPL_H

#include <stdexcept>
#include <cmath>
#include <netcdf>

#include "meanshift_op.h"

namespace m3D { 


    template <typename T>
    void 
    MeanshiftOperation<T>::prime_index(const SearchParameters *params)
    {
        vector<T> x(this->feature_space->rank());
        typename Point<T>::list *sample = this->point_index->search(x,params,NULL);
        delete sample;
    }

    template <typename T>
    vector<T>
    MeanshiftOperation<T>::meanshift(const vector<T> &x,
                                     const SearchParameters *params,
                                     const Kernel<T> *kernel,
                                     const WeightFunction<T> *w,
                                     const bool normalize_shift)
    {
        using namespace utils::vectors;
        
        // Distances are only used in the KNN case.

        vector<T> distances;

        typename Point<T>::list *sample = NULL;

        if ( params->search_type() == SearchTypeRange )
        {
            sample = this->point_index->search( x, params, NULL );

            #if DEBUG_MEANSHIFT_SAMPLING
                RangeSearchParams<T> *rsp = (RangeSearchParams<T> *)params;
                cout << "Sample at " << x << " with search radius " << rsp->bandwidth << " yields " << sample->size() << " points:" << endl;
                for ( size_t i = 0; i < sample->size(); i++ )
                {
                    Point<T> *p = sample->at(i);
                    cout << (i+1) << "\t" << p->values << endl;
                }
            #endif

            #if WRITE_MEANSHIFT_SAMPLES
                static size_t sample_index = 0;
                NetCDFDataStore<T> *ds = (NetCDFDataStore<T> *) this->feature_space->data_store();
                string sample_filename = ds->filename() + "_sample_" + boost::lexical_cast<string>(sample_index++) + ".vtk";
                VisitUtils<T>::write_pointlist_vtk( sample_filename, sample, this->feature_space->coordinate_system->rank() );
                sample_index++;
            #endif

        }
        else
        {
            sample = this->point_index->search( x, params, &distances );
        }

        #if WRITE_MEANSHIFT_WEIGHTS

            static size_t weight_plot_number;

            vector<T> kernel_weights;
            vector<T> variable_weights;
            vector<T> combined_weights;

            bool is_sample_point = false;

            // get the spatial component of x

            vector<T> xs = this->feature_space->spatial_component( x );

            this->feature_space->round_to_resolution( xs, this->feature_space->coordinate_system->resolution() );

            // figure out if the spatial coordinate is part of the
            // vector 'weight_sample_points' in feature-space

            typename vector< vector<T> >::iterator sample_points;

            for ( sample_points = this->feature_space->weight_sample_points.begin();
                 sample_points != this->feature_space->weight_sample_points.end() && !is_sample_point;
                 sample_points++ )
            {
                vector<T> point = *sample_points;

                this->feature_space->round_to_resolution( point, this->feature_space->coordinate_system->resolution() );

                is_sample_point = point == xs;
            }
        #endif

        vector<T> shift( this->feature_space->dimension, 0.0 );

        // If the sample is empty, no shift can be calculated.
        // Returns a shift of 0
        if ( sample->size() == 0 )
        {
            return shift;
        }

        // pre-allocate to avoid constructor/destructor overhead in loop

        vector<T> numerator( this->feature_space->dimension );
        T denominator = 0.0;

        // unroll loop in blocks

        size_t size = sample->size();

        // Range Search

        if ( params->search_type() == SearchTypeRange )
        {
            vector<T> h = ((RangeSearchParams<T> *)  params)->bandwidth;

            if ( w == NULL && kernel != NULL)
            {
                // No weight function

                for ( size_t index = 0; index < size; index++ )
                {
                    T weight = kernel->apply(mahalabonis_distance_sqr(x,sample->at(index)->values,h));
                    denominator += weight;
                    numerator += weight * sample->at(index)->values;

                    #if WRITE_MEANSHIFT_WEIGHTS
                        if ( is_sample_point )
                        {
                            kernel_weights.push_back( weight );
                            variable_weights.push_back( 1.0 );
                            combined_weights.push_back( weight );
                        }
                    #endif

                }
            }
            else if (w != NULL && kernel == NULL)
            {
                // no kernel

                for ( size_t index = 0; index < size; index++ )
                {
                    T weight = w->operator()(sample->at(index));
                    denominator += weight;
                    numerator += weight * sample->at(index)->values;

                    #if WRITE_MEANSHIFT_WEIGHTS
                        if ( is_sample_point )
                        {
                            kernel_weights.push_back( 1.0 );
                            variable_weights.push_back( weight );
                            combined_weights.push_back( weight );
                        }
                    #endif
                }
            }
            else if (kernel == NULL && w == NULL)
            {
                // no weight function AND no kernel

                for ( size_t index = 0; index < size; index++ )
                {
                    denominator += 1.0;
                    numerator += sample->at(index)->values;

                    #if WRITE_MEANSHIFT_WEIGHTS
                        if ( is_sample_point )
                        {
                            kernel_weights.push_back( 1.0 );
                            variable_weights.push_back( 1.0 );
                            combined_weights.push_back( 1.0 );
                        }
                    #endif
                }
            }
            else
            {
                // Both weight function and kernel

                vector<T> mult(x.size(),0.0);

//                #if WITH_OPENMP
//                #pragma omp parallel for schedule(static) // reduction(+:denominator)
//                #endif  
                for ( size_t index = 0; index < size; index++ )
                {
                    T var_weight = w->operator()(sample->at(index));
                    T kernel_weight = kernel->apply(mahalabonis_distance_sqr(x,sample->at(index)->values,h));
                    T weight = kernel_weight * var_weight;

                    denominator += weight;

                    for (int i=0; i<x.size(); i++) 
                        numerator[i] += (weight * sample->at(index)->values[i]);

                    #if WRITE_MEANSHIFT_WEIGHTS
                        if ( is_sample_point )
                        {
                            kernel_weights.push_back( kernel_weight );
                            variable_weights.push_back( var_weight );
                            combined_weights.push_back( weight );
                        }
                    #endif
                }
            }
        }

        // KNN

        else
        {
            throw "Not Implemented";
        }

        // Put it together

#if DEBUG_MEANSHIFT_RESULT_CALCULATION
        cout << endl << "\tnumerator = " << numerator;
        cout << endl << "\tdenominator = " << denominator;

        if ( denominator == 0 )
        {
            cerr << "Denominator is 0!!" << endl;
        }
#endif
        vector<T> dx = numerator / denominator;
        shift = dx - x;

#if DEBUG_MEANSHIFT_RESULT_CALCULATION
        cout << endl << "\tm = " << dx;
        cout << endl << "\tx = " << x;
        cout << endl << "\tshift = " << shift;
        cout << endl;
#endif

#if WRITE_MEANSHIFT_WEIGHTS

        if ( is_sample_point )
        {
            NetCDFDataStore<T> *ds = (NetCDFDataStore<T> *) this->feature_space->data_store();
            string filename = ds->filename() + "_kernel_weights_" + boost::lexical_cast<string>(weight_plot_number) + ".vtk";
            boost::replace_all( filename, "/", "_" );
            boost::replace_all( filename, "..", "" );
            VisitUtils<T>::write_weights( filename, "kernel_weights", sample, kernel_weights );

            filename = ds->filename() + "_variable_weights_" + boost::lexical_cast<string>(weight_plot_number) + ".vtk";
            boost::replace_all( filename, "/", "_" );
            boost::replace_all( filename, "..", "" );
            VisitUtils<T>::write_weights( filename, "variable_weights", sample, variable_weights );

            filename = ds->filename() + "_combined_weights_" + boost::lexical_cast<string>(weight_plot_number) + ".vtk";
            boost::replace_all( filename, "/", "_" );
            boost::replace_all( filename, "..", "" );
            VisitUtils<T>::write_weights( filename, "combined_weights", sample, combined_weights );

            filename = ds->filename() + "_shifts_" + boost::lexical_cast<string>(weight_plot_number) + ".vtk";
            boost::replace_all( filename, "/", "_" );
            boost::replace_all( filename, "..", "" );

            // Cook down to 2D for displaying on 2D plots

            vector< vector<T> > sample_origins;
            vector< vector<T> > sample_shifts;
            vector<T> sample_shift = this->feature_space->spatial_component( shift );
            sample_shifts.push_back( sample_shift );

            xs.push_back(0.0);
            sample_origins.push_back( xs );
            VisitUtils<T>::write_vectors_vtk( filename, sample_origins, sample_shifts );

            weight_plot_number++;
        }
#endif

        delete sample;

        if (normalize_shift)
        {
            shift = this->feature_space->coordinate_system->round_to_grid(shift);
        }

        return shift;
    }
}

#endif
