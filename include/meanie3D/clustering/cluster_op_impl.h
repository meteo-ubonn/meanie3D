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

#ifndef M3D_OPERATION_CLUSTER_IMPL_H
#define M3D_OPERATION_CLUSTER_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/utils.h>

#include <vector>

#include "cluster_task.h"

namespace m3D {

#if WRITE_MODES
    static size_t s_pass_counter = 0;

    template <typename T>
    size_t
    ClusterOperation<T>::pass_counter()
    {
        return s_pass_counter;
    };

    template <typename T>
    void
    ClusterOperation<T>::reset_pass_counter()
    {
        s_pass_counter = 0;
    };

    template <typename T>
    void
    ClusterOperation<T>::increment_pass_counter()
    {
        s_pass_counter++;
    }
#endif

#pragma mark -
#pragma mark Clustering Code

    template <typename T>
    ClusterList<T>
    ClusterOperation<T>::cluster(
            const SearchParameters *params,
            const Kernel<T> *kernel,
            const WeightFunction<T> *weight_function,
            const bool coalesceWithStrongestNeighbour,
            const bool show_progress_bar)
    {
        using namespace m3D::utils::vectors;
        vector<T> resolution;
        if (params->search_type() == SearchTypeRange) {
            RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;
            // Physical grid resolution in the
            // spatial range
            resolution = this->feature_space->coordinate_system->resolution();
            resolution = ((T) 4.0) * resolution;
            // Supplement with bandwidth values for
            // the value range
            for (size_t i = resolution.size(); i < p->bandwidth.size(); i++) {
                resolution.push_back(p->bandwidth[i]);
            }
        } else {
            KNNSearchParams<T> *p = (KNNSearchParams<T> *) params;
            resolution = p->resolution;
        }

#if WRITE_TRAJECTORIES
        static size_t trajectory_number = 0;
#endif  

        if (show_progress_bar) {
            cout << endl << "Creating meanshift vector graph ...";
            start_timer();
            m_progress_bar = new boost::progress_display(this->feature_space->size());
        }

        // Compile list of dimension names

        const CoordinateSystem<T> *cs = this->feature_space->coordinate_system;
        NcFile file(m_data_store->filename(), NcFile::read);
        vector<NcVar> vars(cs->dimension_variables());

        for (size_t i = 0; i < m_data_store->rank(); i++) {
            NcVar variable = file.getVar(m_data_store->variable_names()[i]);
            vars.push_back(variable);
        }
        ClusterList<T> cluster_list(vars, cs->dimensions(),
                m_data_store->filename(),
                m_data_store->get_time_index());

        // Guard against empty feature-space
        if (this->feature_space->points.size() == 0) {
            cout << "Feature space is empty" << endl;
            return cluster_list;
        }

        MeanshiftOperation<T> meanshiftOperator(this->feature_space, this->point_index);
        meanshiftOperator.prime_index(params);

#if WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t index = 0; index < this->feature_space->size(); index++) {
            if (show_progress_bar) {
#if WITH_OPENMP
#pragma omp critical
#endif
                m_progress_bar->operator++();
            }
            // Get the meanshift vector for this point
            typename Point<T>::ptr x = this->feature_space->points[ index ];
            x->shift = meanshiftOperator.meanshift(x->values, params, kernel, weight_function);
            // Extract the spatial component and obtain the grid
            vector<T> spatial_shift = this->feature_space->spatial_component(x->shift);
            x->gridded_shift = this->feature_space->coordinate_system->to_gridpoints(spatial_shift);
        }

        if (show_progress_bar) {
            cout << "done. (" << stop_timer() << "s)" << endl;
            delete m_progress_bar;
            m_progress_bar = NULL;
        }

        // Analyse the graph and create clusters
        cluster_list.aggregate_cluster_graph(this->feature_space, weight_function, coalesceWithStrongestNeighbour, show_progress_bar);

        // Provide fresh ids right away
        m3D::uuid_t uuid = 0;
        ClusterUtils<T>::provideUuids(&cluster_list,uuid);
        m3D::id_t id = 0;
        ClusterUtils<T>::provideIds(&cluster_list,id);

        cout << "Cluster list after aggregation:" << endl;
        cluster_list.print();

        // Replace points with original data ()
        ClusterUtils<T>::replace_points_from_datastore(cluster_list, m_data_store);
        cout << "Cluster list after filtering points:" << endl;
        cluster_list.print();

        // Find margin points (#325)
        ClusterUtils<T>::obtain_margin_flag(cluster_list, this->feature_space);

#if WRITE_BOUNDARIES
        cluster_list.write_boundaries(weight_function, this->feature_space, this->point_index, resolution);
#endif

#if WRITE_MODES
#if WITH_VTK
        size_t min_size = std::numeric_limits<size_t>::max();
        size_t max_size = std::numeric_limits<size_t>::min();
        for (size_t i = 0; i < cluster_list.size(); i++) {
            if (cluster_list[i]->size() < min_size) {
                min_size = cluster_list[i]->size();
            }
            if (cluster_list[i]->size() > max_size) {
                max_size = cluster_list[i]->size();
            }
        }

        NetCDFDataStore<T> *ds = (NetCDFDataStore<T> *) this->feature_space->data_store();
        std::string fn = ds->filename() + "-modes-" + boost::lexical_cast<string>(pass_counter()) + ".vtk";
        VisitUtils<T>::write_cluster_modes_vtk(fn, cluster_list.clusters);

        fn = ds->filename() + "-raw-modes-" + boost::lexical_cast<string>(pass_counter()) + ".vtk";
        VisitUtils<T>::write_modes_vtk(fn, cluster_list.trajectory_endpoints(), cluster_list.trajectory_lengths());
        VisitUtils<T>::write_cluster_modes_vtk(fn, cluster_list.clusters);
#endif
#endif
        return cluster_list;
    }
}

#endif
