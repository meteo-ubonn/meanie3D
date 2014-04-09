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
                       CoordinateSystem<T> *cs,
                       const NcVar &track_variable,
                       Verbosity verbosity)
    {
        if ( verbosity >= VerbosityNormal )
        {
            cout << endl << "Tracking ..." << endl;
            start_timer();
        }


#pragma mark -
#pragma mark Preliminaries        
        
        if ( current->clusters.size() == 0 )
        {
            if ( verbosity >= VerbosityNormal )
                cout << "No current clusters. Skipping tracking." << endl;
            
            return;
        }
        
        // figure out tracking variable index (for histogram)
        
        int index = index_of_first( current->feature_variables, track_variable );
        
        assert( index >= 0 );
        
        size_t tracking_var_index = (size_t) index;
        
        // Check cluster sizes
        
        for (size_t i=0; i < previous->clusters.size(); i++)
        {
            typename Cluster<T>::ptr c = previous->clusters[i];
            
            if (c->points.size()==0)
            {
                cerr << "ERROR: previous cluster " << c->id << " has no points!" << endl;
            }
        }

        for (size_t i=0; i < current->clusters.size(); i++)
        {
            typename Cluster<T>::ptr c = current->clusters[i];
            
            if (c->points.size()==0)
            {
                cerr << "ERROR: current cluster " << c->id << " has no points!" << endl;
            }
        }

        
        // Valid min/max of tracking variable
        
        T valid_min, valid_max;
        
        get_valid_range(track_variable, valid_min, valid_max );
        
        // Find out the highest id from the previous list
        
        size_t highest_id = previous->highest_id;
        
        if (highest_id == NO_ID || (highest_id == 0 && previous->clusters.size() > 0))
        {
            // figure it out from the highest found ID
            
            for (size_t ci=0; ci < previous->clusters.size(); ci++)
            {
                if (previous->clusters[ci]->id > highest_id)
                {
                    highest_id = previous->clusters[ci]->id;
                }
            }
        }
        
        if ( verbosity >= VerbosityDetails )
            cout << "Highest used ID is " << highest_id << endl;
        
        // check time difference and determine displacement restraints
        
        ::units::values::s p_time = ::units::values::s(previous->timestamp);
        ::units::values::s c_time = ::units::values::s(current->timestamp);
        
        // Check c > p
        
        if (p_time >= c_time)
        {
            cerr << "ERROR:previous timestamp not greater than current timestamp!" << endl;
            return;
        }
        
        this->m_deltaT = c_time - p_time;
        
        if (this->m_deltaT > m_max_deltaT)
        {
            cerr << "ERROR:files too far apart in time. Time difference:" << this->m_deltaT << "s, longest accepted:" << m_max_deltaT << "s" << endl;
            return;
        }
        
        current->tracking_time_difference = m_deltaT.get();

        ::units::values::m maxDisplacement = this->m_maxVelocity * this->m_deltaT;
        
        if ( verbosity >= VerbosityDetails )
            printf("max velocity  constraint at %4.1f m/s at deltaT %4.0fs -> dR_max = %7.1fm\n",
                   this->m_maxVelocity.get(), this->m_deltaT.get(), maxDisplacement.get() );

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
        
        ::units::values::meters_per_second overlap_constraint_velocity = m_maxVelocity;
        
        ::units::values::m overlap_constraint_radius = 0.5 * m_deltaT * overlap_constraint_velocity;
        
        // Prepare the newcomers for re-identification
        
        current->erase_identifiers();

        // Get the counts
        
        size_t old_count = previous->clusters.size();
        size_t new_count = current->clusters.size();
        
        if ( verbosity >= VerbosityDetails )
            cout << "Matching " << new_count << " new clusters against " << old_count << " old clusters" << endl;
        
#pragma mark -
#pragma mark Compute correlation matrixes and constraints
        
        typename Matrix<T>::matrix_t rank_correlation = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix< ::units::values::m >::matrix_t midDisplacement = Matrix< ::units::values::m >::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t histDiff = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t sum_prob = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t coverOldByNew = Matrix<T>::create_matrix(new_count,old_count);
        typename Matrix<T>::matrix_t coverNewByOld = Matrix<T>::create_matrix(new_count,old_count);
        
        typename Matrix<T>::flag_matrix_t constraints_satisified = Matrix<T>::create_flag_matrix(new_count, old_count);
        
        int maxHistD = numeric_limits<int>::min();
        ::units::values::m maxMidD = ::units::values::m(numeric_limits<T>::min());
        
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
//                // that a cloud covers 10 pixels in one scan and unit_multiplier in the next. The
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
                    ::units::values::m radius = oldCluster->radius(cs);
                    
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
                
                midDisplacement[n][m] = ::units::values::m(vector_norm(cs->to_meters(dx)));
                
                //
                // Maximum velocity constraint
                //
                
                // TODO: calculate the max displacement in the same dimension
                // as the dimension variables
                
                ::units::values::m displacement = midDisplacement[n][m];
                
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
                
                if (midDisplacement[n][m] > maxMidD)
                {
                    maxMidD = midDisplacement[n][m];
                }
                
                // only if all constraints are passed, the flag is set to true
                
                constraints_satisified[n][m] = true;
            }
        }
        
        // Can't have zeros here
        
        if (maxHistD==0)
        {
            maxHistD = 1;
        }
        
#pragma mark -
#pragma mark Calculate probabilities
        
        if ( verbosity >= VerbosityDetails )
            cout << endl << "-- Correlation Table --" << endl;
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr newCluster = current->clusters[n];
            
            if ( verbosity >= VerbosityDetails )
            {
                cout << endl << "Correlating new Cluster #" << n
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
                    float prob_r = m_dist_weight * erfc( midDisplacement[n][m].get() / maxMidD.get() );
                    
                    float prob_h = m_size_weight * erfc( histDiff[n][m] / maxHistD );
                    
                    float prob_t = m_corr_weight * rank_correlation[n][m];
                    
                    sum_prob[n][m] = prob_t + prob_r + prob_h;
                    
                    if ( verbosity >= VerbosityDetails )
                    {
                        printf("\t<ID#%4lu>:\t(|H|=%5lu)\t\tdR=%4.1f (%5.4f)\t\tdH=%5.4f (%5.4f)\t\ttau=%7.4f (%5.4f)\t\tsum=%6.4f\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                               oldCluster->id,
                               oldCluster->histogram(tracking_var_index,valid_min,valid_max)->sum(),
                               midDisplacement[n][m].get(),
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

#pragma mark -
#pragma mark Matchmaking
        
        ::units::values::meters_per_second velocitySum = ::units::values::meters_per_second(0);
        
        int velocityClusterCount = 0;
        
        float currentMaxProb = numeric_limits<float>::max();
        
        int maxIterations = new_count * old_count;
        
        int iterations = 0;
        
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
            
            if ( new_cluster->id == cfa::NO_ID )
            {
                // new cluster not tagged yet
                
                if ( used_clusters.find(old_cluster) == used_clusters.end() )
                {
                    // old cluster not matched yet

                    // TODO: calculate velocity in the same units as the dimension
                    // variables
                    
                    ::units::values::meters_per_second velocity = midDisplacement[maxN][maxM] / this->m_deltaT;
                
                    velocitySum += velocity;
                        
                    velocityClusterCount++;
                    
                    if ( verbosity >= VerbosityDetails )
                    {
                        printf("pairing new cluster #%4lu / old Cluster ID=#%4lu accepted, velocity %4.1f m/s\n", maxN, old_cluster->id, velocity.get() );
                    }
                        
                    new_cluster->id = old_cluster->id;
                    
                    used_clusters.insert( old_cluster );
                    
                    current->tracked_ids.insert( new_cluster->id );
                }
            }
            
            currentMaxProb = maxProb;
            
            iterations++;
        }

//        // update mean velocity
//        float newMeanVelocity = 0;
//        if (velocityClusterCount>0) newMeanVelocity = velocitySum/velocityClusterCount;
//        if (newMeanVelocity > 2.0) meanVelocity = newMeanVelocity;
//
//        SPLog(@"\ncurrent mean velocity: %4.2f m/s",meanVelocity);
        
#pragma mark -
#pragma mark Tag unmatched clusters

        if ( verbosity >= VerbosityDetails )
            printf("\n-- New ID assignments --\n");
        
        for ( n=0; n < new_count; n++ )
        {
            typename Cluster<T>::ptr c = current->clusters[n];
            
            if ( c->id == cfa::NO_ID )
            {
                c->id = next_id(highest_id);
                
                current->new_ids.insert( c->id );
                
                if ( verbosity >= VerbosityDetails )
                {
                    cout << "#" << n
                         << " at " << c->mode
                         << " with |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum()
                         << " is assigned new ID #" << c->id
                         << endl;
                }
            }
        }

#pragma mark -
#pragma mark Determine drop-outs
        
        bool hadGoners = false;
        
        typename Cluster<T>::list::iterator ci;
        
        if ( verbosity >= VerbosityDetails )
            cout << "\n-- Goners --" << endl;
        
        for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            typename set< typename Cluster<T>::ptr >::iterator f = used_clusters.find( c );
            
            if ( f != used_clusters.end() ) continue;
            
            if ( verbosity >= VerbosityDetails )
                cout << "#" << c->id << " at " << c->mode << " |H|=" << c->histogram(tracking_var_index,valid_min,valid_max)->sum() << endl;
            
            current->dropped_ids.insert(c->id);

            hadGoners = true;
            
        }
        if (!hadGoners)
        {
            cout << "none." << endl;
        }
        
#pragma mark -
#pragma mark Merging
        
        if ( verbosity >= VerbosityDetails )
            cout << endl << "-- Merges --" << endl;
        
        // for each new Cluster check coverage of old Clusters
        
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
                    
                    if ( overlap > this->m_ms_threshold )
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
                
                // collect merged cluster ids

                id_set_t merged_from;
                
                for ( int i=0; i < candidates.size(); i++ )
                    merged_from.insert(previous->clusters[candidates[i]]->id);
                
                if ( verbosity >= VerbosityDetails )
                    cout << "clusters " << merged_from << " seem to have merged into cluster " << new_cluster->id << endl;
                
                // store the given ID
                
                size_t the_id = new_cluster->id;
                
                // check if the ID was new?
                
                id_set_t::const_iterator it = current->new_ids.find(the_id);
                
                if ( it == current->new_ids.end() )
                {
                    if ( maxCover > this->m_msc_threshold )
                    {
                        // if the biggest congruence is at least 75%, continue the track

                        typename Cluster<T>::ptr c = previous->clusters[ largestCandidateIndex ];
                        
                        new_cluster->id = c->id;
                        
                        current->tracked_ids.insert(new_cluster->id);
                        
                        if ( verbosity >= VerbosityDetails )
                            printf("\ttrack ID#%lu continues\n", new_cluster->id );
                    }
                    else
                    {
                        // Otherwise dish out a fresh tag and add to new IDs
                        
                        new_cluster->id = next_id(highest_id);
                        
                        current->new_ids.insert(new_cluster->id);
                        
                        // as well as remove from tracked IDs
                        
                        id_set_t::iterator fi = current->tracked_ids.find(the_id);
                        if (fi!=current->tracked_ids.end())
                        {
                            current->tracked_ids.erase(the_id);
                        }

                        if ( verbosity >= VerbosityDetails )
                            printf("\ttrack ID#%lu ends. new track ID#%lu begins\n", the_id, new_cluster->id );
                    }
                }
                else
                {
                    if ( verbosity >= VerbosityDetails )
                        printf("\ttrack ID#%lu continues\n", new_cluster->id );
                }
                
                // store the merge in the current cluster file
                
                current->merges[new_cluster->id] = merged_from;
            }
        }
        
        if ( ! had_merges )
            if ( verbosity >= VerbosityDetails )
                cout << "none." << endl;
        
#pragma mark -
#pragma mark Splitting
        
        if ( verbosity >= VerbosityDetails )
            cout << endl << "-- Splits --" << endl;
        
        // for each new Cluster check coverage of old Clusters
        
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

                    if ( overlap > this->m_ms_threshold )
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
                
                id_set_t split_into;
                
                for ( int i=0; i < candidates.size(); i++ )
                    split_into.insert(current->clusters[candidates[i]]->id);
                
                if ( verbosity >= VerbosityDetails )
                    cout << "cluster " << old_cluster->id << " seems to have split into clusters " << split_into << endl;
                
                for ( int i=0; i < candidates.size(); i++ )
                {
                    typename Cluster<T>::ptr c = current->clusters[ candidates[i] ];
                    
                    // check if the tracked id's contains c->id
                    
                    cfa::id_t the_id = c->id;
                    
                    id_set_t::iterator it = current->tracked_ids.find(the_id);

                    // If the largest one of the new clusters is at least 75% of the
                    // size of new previous cluster, the ID of the previous cluster
                    // is continued in the largest candidate
                    
                    bool continueID = (candidates[i] == largestCandidateIndex) && (maxCover > this->m_msc_threshold);
                    
                    if (continueID)
                    {
                        c->id = old_cluster->id;
                    
                        current->tracked_ids.insert(c->id);
                        
                        if ( verbosity >= VerbosityDetails )
                            printf("\ttrack ID#%lu continues\n", c->id );
                    }
                    else
                    {
                        if ( it != current->tracked_ids.end() )
                        {
                            // remove from tracked id's, re-tag and add to new IDs
                            
                            current->tracked_ids.erase(it);
                            
                            c->id = next_id(highest_id);
                            
                            if ( verbosity >= VerbosityDetails )
                                printf("\tnew track ID#%lu begins\n", c->id );
                            
                            current->new_ids.insert( c->id );
                        }
                    }
                }
                
                current->splits[old_cluster->id] = split_into;
            }
        }
        
        if ( ! had_splits)
            if ( verbosity >= VerbosityDetails )
                cout << "none." << endl;

#pragma mark -
#pragma mark Store tracking data
        
        if ( verbosity >= VerbosityNormal )
        {
            cout << "done. (Correlating " << (new_count * old_count) << " objects took " << stop_timer() << " seconds)" << endl;
        }
        
        current->tracking_performed = true;
        
        current->highest_id = highest_id;
    }
    
}

#endif
