#ifndef _M3D_ClusterUtils_Impl_H_
#define _M3D_ClusterUtils_Impl_H_

#include <cf-algorithms/cf-algorithms.h>

#include <meanie3D/utils.h>

namespace m3D {

    template <typename T>
    ClusterUtils<T>::ClusterUtils(float merge_threshold)
    : m_merge_threshold(merge_threshold) {}

    template <typename T>
    void
    ClusterUtils<T>::filter_with_previous_clusters(typename ClusterList<T>::ptr previous,
                                                   typename ClusterList<T>::ptr current,
                                                   WeightFunction<T> *weight_function,
                                                   const Verbosity verbosity)
    {
        if ( verbosity >= VerbosityNormal )
            cout << endl << "-- Filtering with previous results --" << endl;
        
        // sanity check
        
        if ( current->clusters.size() == 0 || previous->clusters.size() == 0)
        {
            if ( verbosity >= VerbosityNormal )
            {
                cout << "Nothing to do." << endl;
                return;
            }
        }
        
        size_t new_count = current->clusters.size();
        size_t old_count = previous->clusters.size();
        ::cfa::meanshift::id_t current_id = current->clusters[current->clusters.size()-1]->id;
        
        typename Matrix<T>::matrix_t coverOldByNew = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t coverNewByOld = Matrix<T>::create_matrix(new_count,old_count);
        
        // compute overlap
        
        size_t n,m;
        
        // check out the values for midDisplacement, histSizeDifference and kendall's tau
        
        for ( m=0; m < old_count; m++ )
        {
            typename Cluster<T>::ptr oldCluster = previous->clusters[m];

            for ( n=0; n < new_count; n++ )
            {
                typename Cluster<T>::ptr newCluster = current->clusters[n];

                // Calculate overlap
                
                T overlap = newCluster->percent_covered_by( oldCluster );
                
                if (overlap > 0)
                {
                    printf("old #%4llu with new #%4llu overlap = %3.2f\n", oldCluster->id, newCluster->id, overlap);
                }
                
                coverNewByOld[n][m] = overlap;
            }
        }
        
        // for each new blob check coverage of old blobs
        
        typedef set< typename Cluster<T>::ptr > cset_t;
        
        cset_t erased, merged;
        
        for ( m=0; m < previous->clusters.size(); m++ )
        {
            typename Cluster<T>::ptr old_cluster = previous->clusters[m];
            
            vector<size_t> candidates;
            
            // figure out the largest candidate
            // TODO: do we still need this here?
            
            for ( n=0; n < current->clusters.size(); n++ )
            {
                T overlap = coverNewByOld[n][m];
                
                if ( overlap > 0.33 )
                {
                    candidates.push_back(n);
                }
            }
            
            if (candidates.size() > 1 )
            {
                if ( verbosity >= VerbosityNormal )
                {
                    printf("Old cluster ID#%4llu seems to have split into new clusters IDs ", old_cluster->id);
                    for ( int i=0; i < candidates.size(); i++ )
                    {
                        typename Cluster<T>::ptr c = current->clusters[candidates[i]];
                        printf("#%llu ",c->id );
                    }
                    cout << endl;
                }
                
                typename Cluster<T>::ptr merged_cluster = new Cluster<T>();
                
                merged_cluster->id = ++current_id;
                
                vector<T> mode(previous->feature_variables.size(),0.0);
                
                // merge with those candidates, that are direct neighbours
                // neglect those, that have no direct boundary with any of
                // the other candidates
                
                size_t num_picked_candidates = 0;
                
                for ( int i = 0; i < candidates.size(); i++ )
                {
                    bool have_boundary = true;
                    
//                    typename Cluster<T>::ptr ci = current->clusters[candidates[i]];
//                    
//                    for ( int j = 0; j < candidates.size() && !have_boundary; j++ )
//                    {
//                        if ( i==j ) continue;
//                        
//                        typename Cluster<T>::ptr cj = current->clusters[candidates[j]];
//                        
//                        have_boundary = current->are_neighbours(ci,cj);
//                    }
                    
                    if (have_boundary)
                    {
                        typename Cluster<T>::ptr c = current->clusters[ candidates[i] ];
                        
                        merged_cluster->add_points(c->points);
                        
                        mode += c->mode;
                        
                        erased.insert(c);
                        
                        num_picked_candidates++;
                    }
                }
                
                if (num_picked_candidates > 1)
                {
                    merged_cluster->mode = mode / ((T)num_picked_candidates);
                
                    merged.insert(merged_cluster);
                }
                else
                {
                    delete merged_cluster;
                }
            }
        }
        
        // remove the erased ones and add the merged ones
        
        typename cset_t::iterator ci;
        
        for (ci=erased.begin(); ci!=erased.end(); ci++)
        {
            typename Cluster<T>::ptr c = *ci;
            
            current->clusters.erase(find(current->clusters.begin(), current->clusters.end(), c));
            
            delete c;
        }
        
        for (ci=merged.begin(); ci!=merged.end(); ci++)
        {
            typename Cluster<T>::ptr c = *ci;
            
            current->clusters.push_back(c);
        }
    }
}

#endif
