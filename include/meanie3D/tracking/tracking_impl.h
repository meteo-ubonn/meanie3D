#ifndef _M3D_TRACKING_CLASS_Impl_H_
#define _M3D_TRACKING_CLASS_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>
#include <set>

#include <numericalrecipes/nr.h>

#include <meanie3D/types/cluster.h>
#include <meanie3D/tracking/tracking.h>
#include <meanie3D/utils.h>

namespace m3D {
    
    using namespace utils;
    using namespace netCDF;

    template <typename T>
    void
    Tracking<T>::track(typename ClusterList<T>::ptr previous,
                       typename ClusterList<T>::ptr current,
                       const NcVar &track_variable,
                       Verbosity verbosity)
    {
        if ( verbosity >= VerbosityNormal )
            cout << endl << "-- Tracking --" << endl;
        
        // sanity check
        
        if ( current->clusters.size() == 0 )
        {
            if ( verbosity >= VerbosityNormal )
            {
                cout << "No current clusters. Skipping tracking." << endl;
                
                return;
            }
        }
        
        // figure out tracking variable index
        
        int index = index_of_first( current->feature_variables, track_variable );
        
        assert( index >= 0 );
        
        size_t tracking_var_index = (size_t) index;
        
        
        // Valid min/max
        
        T valid_min, valid_max;
        
        track_variable.getAtt("valid_min").getValues( &valid_min );

        track_variable.getAtt("valid_max").getValues( &valid_max );
        
        // Find the highest id
        
        size_t current_id = 0;

        // Find out the highest id from the previous list
        
        typename Cluster<T>::list::iterator ci;
        
        for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            if ( c->id > current_id )
            {
                current_id = c->id;
            }
        }
        
        current_id++;
        
        if ( verbosity >= VerbosityNormal )
            cout << "next available id is " << current_id << endl;
        
        // check time difference and determine displacement restraints
        
        timestamp_t p_time = previous->timestamp;
        
        timestamp_t c_time = current->timestamp;
        
        // Check c > p
        
        if (p_time >= c_time)
        {
            cerr << "ERROR:previous timestamp not greater than current timestamp!" << endl;
            return;
        }
        
        this->m_deltaT = boost::numeric_cast<double>(c_time) - boost::numeric_cast<double>(p_time);
        
        if (this->m_deltaT > m_max_deltaT)
        {
            cerr << "ERROR:files too far apart in time. Time difference:" << this->m_deltaT << "s, longest accepted:" << m_max_deltaT << "s" << endl;
            return;
        }
        

        T maxDisplacement = this->m_maxVelocity * this->m_deltaT;
        
        if ( verbosity >= VerbosityNormal )
            printf("max velocity  constraint at %4.1f m/s at deltaT %4.0fs -> dR_max = %7.1fm\n",
                   this->m_maxVelocity, this->m_deltaT, maxDisplacement );

        // TODO: mean velocity constraint
        
//        T maxMeanVelocityDisplacement = 0.0;
//        
//        if (m_useMeanVelocityConstraint)
//        {
//            m_useMeanVelocityConstraint = deltaT * meanVelocity * m_meanVelocitySecurityPercentage;
//            
//            printf("mean velocity constraint at %4.2f m/s, dR_max = %7.1f",
//                   meanVelocity * m_meanVelocitySecurityPercentage, maxMeanVelocityDisplacement );
//        }
        
//        SPLog(@"\n");
        
        // Minimum object radius for overlap constraint
        // TODO: this point needs to be replaced with velocity
        // vectors from SMV/RMV estimates
        
        T overlap_constraint_velocity = m_maxVelocity;
        
        T overlap_constraint_radius = 0.5 * m_deltaT * overlap_constraint_velocity;
        
        // Prepare the newcomers for re-identification
        
        current->erase_identifiers();

        // Get the counts
        
        size_t old_count = previous->clusters.size();
        size_t new_count = current->clusters.size();
        
        if ( verbosity >= VerbosityNormal )
            cout << "Matching " << new_count << " new clusters against " << old_count << " old clusters" << endl;
        
        // Create result matrices
        
        typename Matrix<T>::matrix_t rank_correlation = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t midDisplacement = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t histDiff = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t sum_prob = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t coverOldByNew = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t coverNewByOld = Matrix<T>::create_matrix(new_count,old_count);
        
        typename Matrix<T>::flag_matrix_t constraints_satisified = Matrix<T>::create_flag_matrix(new_count, old_count);
        
        T maxHistD = numeric_limits<T>::min();
        
        T maxMidD = numeric_limits<T>::min();
        
        // compute correlation matrix
        
        int n,m;
        
        // check out the values for midDisplacement, histSizeDifference and kendall's tau
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current->clusters[n];
            
            typename Histogram<T>::ptr newHistogram = newCluster->histogram(tracking_var_index,valid_min,valid_max);
            
            for ( m=0; m < old_count; m++ )
            {
                // set all constraint flags to false to start with
                
                constraints_satisified[n][m] = false;
                
                // Displacement
                
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];
                
                // Histogram Size
                
                typename Histogram<T>::ptr oldHistogram = oldCluster->histogram(tracking_var_index,valid_min,valid_max);
//
//                size_t max_size = max( newHistogram->sum(), oldHistogram->sum() );
//                
//                histDiff[n][m] = (max_size==0) ? 1.0 : (abs( (T)newHistogram->sum() - (T)oldHistogram->sum() ) / ((T)max_size) );
//                
//                //
//                // Size deviation overlap constraint
//                //
//                
//                // Processes in nature develop within certain bounds. It is not possible
//                // that a cloud covers 10 pixels in one scan and 1000 in the next. The
//                // size deviation constraint is created to prohibit matches between objects,
//                // which vary too much in size
//                
//                T max_H = (T) max(oldHistogram->sum(), newHistogram->sum());
//                
//                T min_H = (T) min(oldHistogram->sum(), newHistogram->sum());
//                
//                T size_deviation = max_H / min_H;
//                
//                if ( size_deviation > m_max_size_deviation ) continue;
                
                //
                // Overlap Constraint
                //
                
                // if the object is so big, that overlap is required at the given advection velocity
                // then check if that is the case. If no overlap exists, prohibit the match by setting
                // the constraint to false. If no overlap is required, the constraint is simply set to
                // true, thus allowing a match.
                
                // Note: radius is calculated in kilometres.
                // TODO: calculate overlap constraint radius in the same dimension
                // as the dimension variables!
                
                coverOldByNew[n][m] = oldCluster->percent_covered_by( newCluster );

                coverNewByOld[n][m] = newCluster->percent_covered_by( oldCluster );

                if (m_useOverlapConstraint)
                {
                    T radius = 1000 * oldCluster->radius();
                    
                    bool requires_overlap = (radius >= overlap_constraint_radius);
                    
                    
                    bool overlap_constraint_satisfied = true;
                    
                    if (requires_overlap)
                    {
                        overlap_constraint_satisfied = coverOldByNew[n][m] > 0.0;
                    }
                    
                    if ( !overlap_constraint_satisfied ) continue;
                }
                
                // calculate average mid displacement
                
                vector<T> oldCenter = oldCluster->geometrical_center(current->dimensions.size());
                vector<T> newCenter = newCluster->geometrical_center(current->dimensions.size());
                
                vector<T> dx = newCenter - oldCenter;
                
                midDisplacement[n][m] = vector_norm(dx);
                
                //
                // Maximum velocity constraint
                //
                
                // TODO: calculate the max displacement in the same dimension
                // as the dimension variables
                
                T displacement = 1000.0 * midDisplacement[n][m];
                
                if ( displacement > maxDisplacement ) continue;


                // Histogram Correlation
                
                // calculate spearman's tau
                
                if ( m_corr_weight != 0.0 )
                {
                    rank_correlation[n][m] = newHistogram->correlate_kendall( oldHistogram );
                }

                // track maxHistD and maxMidD
                
                if ( histDiff[n][m] > maxHistD )
                {
                    maxHistD = histDiff[n][m];
                }
                
                if (midDisplacement[n][m]>maxMidD)
                {
                    maxMidD = midDisplacement[n][m];
                }
                
                // only if all constraints are passed, the flag is set to true
                
                constraints_satisified[n][m] = true;
            }
        }
        
        // Can't have zeros here
        
        if (maxHistD==0) maxHistD = 1;
        
        // normalise maxMidD with maximum possible distance
        // why?
        
        // maxMidD = sqrt( 2*200*200 );

        // calculate the final probability
        
        if ( verbosity >= VerbosityDetails )
            cout << endl << "-- Correlation Table --" << endl;
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current->clusters[n];
            
            if ( verbosity >= VerbosityNormal )
            {
                cout << endl << "Correlating new Blob #" << n
                << " (histogram size " << newCluster->histogram(tracking_var_index,valid_min,valid_max)->sum() << ")"
                << " at " << newCluster->mode << endl;
            }

            for ( m=0; m < old_count; m++ )
            {
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];
                
                // Only calculate values for pairs, that satisfy the
                // overlap constraint
                
                if ( constraints_satisified[n][m] )
                {
                    float prob_r = m_dist_weight * erfc( midDisplacement[n][m] / maxMidD );
                    float prob_h = m_size_weight * erfc( histDiff[n][m] / maxHistD );
                    float prob_t = m_corr_weight * rank_correlation[n][m];
                    sum_prob[n][m] = prob_t + prob_r + prob_h;
                    
                    if ( verbosity >= VerbosityNormal )
                    {
                        printf("\t<ID#%4lu>:\t(|H|=%5lu)\t\tdR=%4.1f (%5.4f)\t\tdH=%5.4f (%5.4f)\t\ttau=%7.4f (%5.4f)\t\tsum=%6.4f\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                               oldCluster->id,
                               oldCluster->histogram(tracking_var_index,valid_min,valid_max)->sum(),
                               midDisplacement[n][m],
                               prob_r,
                               histDiff[n][m],
                               prob_h,
                               rank_correlation[n][m],
                               prob_t,
                               sum_prob[n][m],
                               coverOldByNew[n][m],
                               coverNewByOld[n][m]);
                    }
                }
                else
                {
                    if ( verbosity >= VerbosityDetails )
                        printf("\t<ID#%4lu>:\toverlap, size or max velocity constraints violated\n", oldCluster->id);
                }
            }
        }
        
        float velocitySum = 0;
        
        int velocityBlobCount = 0;
        
        float currentMaxProb = numeric_limits<float>::max();
        
        int maxIterations = new_count * old_count;
        
        int iterations = 0;
        
        // NSMutableSet* usedOldBlobs = [[NSMutableSet alloc] init];
        
        set< typename Cluster<T>::ptr > used_clusters;
        
        if ( verbosity >= VerbosityDetails )
            printf("\n-- Matching Results --\n");
        
        while ( iterations <= maxIterations )
        {
            // find the highest significant correlation match
            
            float maxProb = numeric_limits<float>::min();
            
            size_t maxN=0, maxM=0;
            
            bool matchFound = false;
            
            for ( n=0; n < new_count; n++ )
            {
                for ( m=0; m < old_count; m++ )
                {
                    // only consider pairings, where the overlap constraint is satisfied
                    
                    if ( constraints_satisified[n][m] && (sum_prob[n][m] > maxProb) && (sum_prob[n][m] < currentMaxProb) )
                    {
                        maxProb = sum_prob[n][m];

                        maxN = n;
                        
                        maxM = m;
                        
                        matchFound = true;
                    }
                }
            }
            
            if (!matchFound) break;
            
            // take pick, remove paired candidates from arrays
            
            typename Cluster<T>::ptr new_cluster = current->clusters[maxN];

            typename Cluster<T>::ptr old_cluster = previous->clusters[maxM];
            
            if ( new_cluster->id == Cluster<T>::NO_ID )
            {
                // new cluster not tagged yet
                
                if ( used_clusters.find(old_cluster) == used_clusters.end() )
                {
                    // old cluster not matched yet

                    // TODO: calculate velocity in the same units as the dimension
                    // variables
                    
                    float velocity = 1000 * midDisplacement[maxN][maxM] / this->m_deltaT;
                
                    velocitySum += velocity;
                        
                    velocityBlobCount++;
                    
                    if ( verbosity >= VerbosityNormal )
                    {
                        printf("pairing new blob #%4lu / old blob ID=#%4lu accepted, velocity %4.1f m/s\n", maxN, old_cluster->id, velocity );
                    }
                        
                    new_cluster->id = old_cluster->id;
                    
                    used_clusters.insert( old_cluster );
                    
                    current->tracked_ids.push_back( new_cluster->id );
                }
            }
            else
            {
                if ( verbosity >= VerbosityNormal )
                {
                    printf("pairing new blob #%4lu / old blob ID=#%4lud, cluster was already tagged as #%4lu",
                           maxN, old_cluster->id, new_cluster->id );
                }

            }
            
            currentMaxProb = maxProb;
            
            iterations++;
        }

//        // update mean velocity
//        float newMeanVelocity = 0;
//        if (velocityBlobCount>0) newMeanVelocity = velocitySum/velocityBlobCount;
//        if (newMeanVelocity > 2.0) meanVelocity = newMeanVelocity;
//
//        SPLog(@"\ncurrent mean velocity: %4.2f m/s",meanVelocity);
        
        // all new blobs without id now get one
        
        if ( verbosity >= VerbosityNormal )
            printf("\n-- New ID assignments --\n");
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr c = current->clusters[n];
            
            if ( c->id == Cluster<T>::NO_ID )
            {
                c->id = current_id;
                
                current->new_ids.push_back( c->id );
                
                if ( verbosity >= VerbosityNormal )
                {
                    cout << "#" << n
                         << " at " << c->mode
                         << " with |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum()
                         << " is assigned new ID #" << c->id
                         << endl;
                }
                
                current_id++;
            }
        }

        // list those old blobs not being lined up again
        
        bool hadGoners = false;
        
        if ( verbosity >= VerbosityNormal )
            cout << "\n-- Goners --" << endl;
        
        for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            typename set< typename Cluster<T>::ptr >::iterator f = used_clusters.find( c );
            
            if ( f != used_clusters.end() ) continue;
            
            if ( verbosity >= VerbosityNormal )
                cout << "#" << c->id << " at " << c->mode << " |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum() << endl;
            
            current->dropped_ids.push_back(c->id);

            hadGoners = true;
            
        }
        if (!hadGoners)
        {
            cout << "none." << endl;
        }
        
        // handle merging
        
        if ( verbosity >= VerbosityNormal )
            cout << endl << "-- Merges --" << endl;
        
        // for each new blob check coverage of old blobs
        
        bool had_merges = false;
        
        for ( n=0; n < current->clusters.size(); n++ )
        {
            vector<float> candidates;
            
            typename Cluster<T>::ptr new_cluster = current->clusters[n];
            
            // keep track of the largest candidate, that is: the cluster from
            // the previous series, that covers most of the area of the current
            // cluster
            
            T maxCover = 0.0;
            
            size_t largestCandidateIndex = 0;

            for ( m=0; m < previous->clusters.size(); m++ )
            {
                if (constraints_satisified[n][m])
                {
                    T overlap = coverOldByNew[n][m];
                    
                    if ( overlap > this->m_merge_threshold )
                    {
                        candidates.push_back(m);
                        
                        if (coverNewByOld[n][m] > maxCover)
                        {
                            maxCover = coverNewByOld[n][m];
                            
                            largestCandidateIndex = m;
                        }
                    }
                }
            }
            
            if (candidates.size() > 1 )
            {
                had_merges = true;

                // remove the merged old blobs from the matching process
                // and print them out.
                
                if ( verbosity >= VerbosityNormal )
                {
                    cout << "Blobs ";
                
                    for ( int i=0; i < candidates.size(); i++ )
                    {
                        typename Cluster<T>::ptr c = previous->clusters[ candidates[i] ];
                        
                        printf("#%lu ",c->id);
                    }
                    
                    printf("seem to have merged into blob #%lu\n", new_cluster->id );
                }
                
                // store the given ID
                
                size_t the_id = new_cluster->id;
                
                // check if the ID was new?
                
                int index = index_of_first<id_t>( current->new_ids, the_id );
                
                if ( index < 0 )
                {
                    // not new => was tracked
                    
                    // if the biggest congruence is at least 75%, continue
                    // the track
                    
                    if ( maxCover > 0.75 )
                    {
                        typename Cluster<T>::ptr c = previous->clusters[ largestCandidateIndex ];
                        
                        new_cluster->id = c->id;
                        
                        if ( verbosity >= VerbosityNormal )
                            printf("Re-assigned ID#%lu\n", new_cluster->id );
                    }
                    else
                    {
                        new_cluster->id = ++current_id;
                        
                        if ( verbosity >= VerbosityNormal )
                            printf("Assigned new ID#%lu\n", new_cluster->id );
                        
                    }
                    
                    // tag fresh and add to new IDs

                    current->new_ids.push_back(new_cluster->id);
                    
                    // and remove from tracked IDs

                    size_t index = index_of_first<id_t>( current->tracked_ids, the_id );
                    
                    current->tracked_ids.erase( current->tracked_ids.begin() + index );
                }
            }
        }
        
        if ( ! had_merges )
            if ( verbosity >= VerbosityNormal )
                cout << "none." << endl;
        
        // handle merging
        
        if ( verbosity >= VerbosityNormal )
            cout << endl << "-- Splits --" << endl;
        
        // for each new blob check coverage of old blobs
        
        bool had_splits = false;
        
        for ( m=0; m < previous->clusters.size(); m++ )
        {
            typename Cluster<T>::ptr old_cluster = previous->clusters[m];
            
            vector<float> candidates;
            
            T maxCover = 0.0;

            size_t largestCandidateIndex = 0;

            for ( n=0; n < current->clusters.size(); n++ )
            {
                if (constraints_satisified[n][m])
                {
                    T overlap = coverNewByOld[n][m];

                    if ( overlap > this->m_merge_threshold )
                    {
                        candidates.push_back(n);
                        
                        if (coverOldByNew[n][m] > maxCover)
                        {
                            maxCover = coverOldByNew[n][m];
                            
                            largestCandidateIndex = n;
                        }
                    }
                }
            }

            if (candidates.size() > 1 )
            {
                had_splits = true;
                
                if ( verbosity >= VerbosityNormal )
                    printf("Blob ID#%lu seems to have split into blobs ", old_cluster->id);
                
                for ( int i=0; i < candidates.size(); i++ )
                {
                    typename Cluster<T>::ptr c = current->clusters[ candidates[i] ];
                    
                    // check if the new blob's IDs are in tracked IDs and
                    // retag / remove them
                    
                    printf("#%lu",c->id );

                    // check if the tracked id's contains c->id
                    
                    int index = index_of_first<id_t>( current->tracked_ids, (size_t)c->id );

                    // If the largest one of the new clusters is at least 75% of the
                    // size of new previous cluster, the ID of the previous cluster
                    // is continued in the largest candidate
                    
                    bool continueID = (candidates[i] == largestCandidateIndex) && (maxCover > 0.75);
                    
                    if (continueID)
                    {
                        c->id = old_cluster->id;
                    
                        if (index < 0)
                        {
                            current->tracked_ids.push_back(c->id);
                        }
                    }
                    else
                    {
                        if ( index >= 0 )
                        {
                            // remove from tracked id's, re-tag and add to new IDs
                            
                            current->tracked_ids.erase( current->tracked_ids.begin() + index );
                            
                            c->id = current_id++;
                            
                            if ( verbosity >= VerbosityNormal )
                                printf(" (re-tagged as #%lu)",c->id);
                            
                            current->new_ids.push_back( c->id );
                        }
                    }
                    
                    if (i<candidates.size()-1)
                        cout << ",";
                }
            
                cout << endl;
            }
        }
        
        if ( ! had_splits)
            if ( verbosity >= VerbosityNormal )
                cout << "none." << endl;
        
        current->tracking_performed = true;
    }
    
}

#endif
