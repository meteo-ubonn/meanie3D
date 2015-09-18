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

#ifndef M3D_TRACKING_CLASS_Impl_H
#define M3D_TRACKING_CLASS_Impl_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering.h>
#include <meanie3D/numericalrecipes.h>
#include <meanie3D/utils.h>

#include <vector>
#include <netcdf>
#include <set>
#include <utility>
#include <algorithm>

#include "tracking.h"

namespace m3D {
    
    template <typename T>
    void
    Tracking<T>::advect_clusters(typename ClusterList<T>::ptr clusters) {
        using namespace utils::vectors;
        
        for (size_t i=0; i < clusters->size(); i++) {
            typename Cluster<T>::ptr c = clusters->operator [](i);
            if (vector_norm<T>(c->displacement) != 0) {
                // shift all points by the displacement vector
                for (size_t i=0; i<c->size(); i++) {
                    typename Point<T>::ptr p = c->get_points()[i];
                    // shift coordinate
                    p->coordinate += c->displacement;
                    // TODO: shift gridded coordinate
                    for (size_t ci=0; ci < c->displacement.size(); ci++) {
                        // shift feature space point spatial range
                       p->values[ci]+= c->displacement[ci];
                    }
                }
            }
        }
    }

    template <typename T>
    void
    Tracking<T>::track(typename ClusterList<T>::ptr previous,
                       typename ClusterList<T>::ptr current,
                       const CoordinateSystem<T> *cs,
                       const std::string *tracking_variable_name,
                       Verbosity verbosity)
    {
        using utils::Matrix;
        using namespace utils;
        using namespace utils::vectors;
        
#pragma mark -
#pragma mark Preliminaries        
        
        bool skip_tracking = false;

        if ( current->clusters.size() == 0 ) {
            if ( verbosity >= VerbosityNormal ) {
                cout << endl << "No current clusters. Skipping tracking." << endl;
            }
            skip_tracking = true;
        }

        // Find out the highest id from the previous list

        m3D::id_t highest_id = previous->highest_id;

        if (highest_id == NO_ID || (highest_id == 0 && previous->clusters.size() > 0)) {
            // figure it out from the highest found ID
            highest_id = 0;
            for (size_t ci=0; ci < previous->clusters.size(); ci++) {
                if (previous->clusters[ci]->id > highest_id) {
                    highest_id = previous->clusters[ci]->id;
                }
            }
        }
        
        if ( verbosity >= VerbosityDetails ) {
            cout << "Highest used ID is " << highest_id << endl;
        }
        
        if (!skip_tracking) {
            // Check cluster sizes
            
            for (size_t i=0; i < previous->clusters.size(); i++) {
                typename Cluster<T>::ptr c = previous->clusters[i];
                if (c->size()==0) {
                    cerr << "ERROR: previous cluster " << c->id << " has no points!" << endl;
                }
            }
            
            for (size_t i=0; i < current->clusters.size(); i++) {
                typename Cluster<T>::ptr c = current->clusters[i];
                if (c->size()==0) {
                    cerr << "ERROR: current cluster " << c->id << " has no points!" << endl;
                }
            }

            // figure out tracking variable index (for histogram)
            int tracking_var_index = -1;
            std::string tracking_variable;
            T valid_min, valid_max;
            if (tracking_variable_name != NULL) {
                for (size_t i=0; i<current->feature_variables.size(); i++) {
                    if (current->feature_variables[i].getName() == *tracking_variable_name) {
                        tracking_variable = current->feature_variables[i].getName();
                        tracking_var_index = i;
                        break;
                    }
                }
                
                if (tracking_var_index < 0) {
                    cerr << "FATAL:Illegal variable name " << *tracking_variable_name << " for histogram comparison " << endl;
                    exit(EXIT_FAILURE);
                }
                
                // Valid min/max of tracking variable
                utils::netcdf::get_valid_range(current->ncFile->getVar(tracking_variable), valid_min, valid_max );
            }

            // check time difference and determine displacement restraints

            ::units::values::s p_time = ::units::values::s(previous->timestamp);
            ::units::values::s c_time = ::units::values::s(current->timestamp);

            // Check c > p

            if (p_time >= c_time) {
                cerr << "ERROR:previous timestamp not greater than current timestamp!" << endl;
                return;
            }

            this->m_deltaT = c_time - p_time;

            if (this->m_deltaT > m_max_deltaT) {
                cerr << "ERROR:files too far apart in time. Time difference:"
                     << this->m_deltaT << "s, longest accepted:"
                     << m_max_deltaT << "s" << endl;
                return;
            }
            
            // set precision for outputting correlation tables
            std::streamsize prec = std::cout.precision();
            std::cout.precision(2);
            cout << setiosflags(ios::fixed);
            
            // calculate max displacement based on time difference
            // and max velocity
            
            current->tracking_time_difference = m_deltaT.get();
            ::units::values::m maxDisplacement = this->m_maxVelocity * this->m_deltaT;

            if ( verbosity >= VerbosityDetails ) {
                cout << "Max velocity constraint at " 
                     << this->m_maxVelocity.get() << "m/s"
                     << " and dT=" << this->m_deltaT.get()
                     << " results in dRmax=" << maxDisplacement.get() 
                    << endl;
            }

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

            // Minimum object radius for overlap constraint
            // TODO: this point needs to be replaced with velocity
            // vectors from SMV/RMV estimates

            ::units::values::meters_per_second overlap_constraint_velocity = m_maxVelocity;
            ::units::values::m overlap_constraint_radius = 0.5 * m_deltaT * overlap_constraint_velocity;

            // Advect the previous clusters by their displacement 
            // vectors to make the results more accurate
            if (this->m_useDisplacementVectors) {
                this->advect_clusters(previous);
            }
            
            // Prepare the newcomers for re-identification
            current->erase_identifiers();

            // Get the counts

            size_t old_count = previous->clusters.size();
            size_t new_count = current->clusters.size();
            
            if ( verbosity >= VerbosityNormal ) {
                cout << "Calculating preliminary data  ..." << flush;
                start_timer();
            }

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

            // Index the cluster lists for quick overlap calculations
            
            // need to label the clusters with a unique id before we
            // index them, otherwise the area calculations will be nonsense
            for (size_t n=0; n < current->clusters.size(); n++) {
                current->clusters[n]->id = n;
            }
            ClusterIndex<T> new_index(current->clusters, cs->get_dimension_sizes());
            ClusterIndex<T> old_index(previous->clusters, cs->get_dimension_sizes());

            // Employ a mapping to parallelize
            
            vector<size_t> index_dims;
            index_dims.push_back(new_count);
            index_dims.push_back(old_count);
            LinearIndexMapping mapping(index_dims);
            
//            #if WITH_OPENMP
//            #pragma omp parallel for schedule(dynamic)
//            #endif
            for (size_t idx = 0; idx < mapping.size(); idx++) {
                vector<int> index_pair = mapping.linear_to_grid(idx);
                
                int n = index_pair[0];
                int m = index_pair[1];

                // set all constraint flags to false to start with
                constraints_satisified[n][m] = false;

                typename Cluster<T>::ptr newCluster = current->clusters[n];
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];

                //
                // Overlap Constraint
                //

                // if the object is so big, that overlap is required at the given advection velocity
                // then check if that is the case. If no overlap exists, prohibit the match by setting
                // the constraint to false. If no overlap is required, the constraint is simply set to
                // true, thus allowing a match.

                // Note: radius is calculated in kilometers.
                // TODO: calculate overlap constraint radius in the same dimension
                // as the dimension variables!

                bool overlap_constraint_satisfied = true;

//                #if WITH_OPENMP
//                #pragma omp critical
//                #endif
                {
                    // How many percent of old cluster's points are shared by
                    // the new cluster?
                    T old_covered_by_new = new_index.occupation_ratio(oldCluster,newCluster);
                    coverOldByNew[n][m] = old_covered_by_new;
                    
                    // How many percent of the new cluster's points are shared
                    // by the old cluster
                    T new_covered_by_old = old_index.occupation_ratio(newCluster,oldCluster);
                    coverNewByOld[n][m] = new_covered_by_old;
                                                        
                    if (m_useOverlapConstraint) {
                        ::units::values::m radius = oldCluster->radius(cs);
                        bool requires_overlap = (radius >= overlap_constraint_radius);
                        if (requires_overlap) {
                            overlap_constraint_satisfied = coverOldByNew[n][m] > 0.0;
                        }
                    }
                    
                    if (overlap_constraint_satisfied) {
                        // calculate average mid displacement, that is the 
                        // distance between the two clusters geometrical
                        // centers. 

                        vector<T> oldCenter,newCenter;
                        oldCenter = oldCluster->geometrical_center();
                        newCenter = newCluster->geometrical_center();

                        vector<T> dx = newCenter - oldCenter;
                        midDisplacement[n][m] = ::units::values::m(vector_norm(cs->to_meters(dx)));
                    }
                }
                
                if (!overlap_constraint_satisfied && verbosity > VerbosityNormal) {
                    cout << "Pairing old #" << oldCluster->id 
                            << " with new cluster #" << n 
                            << " rejected because of violation of overlap constraint"
                            << endl;
                    continue;
                }
                
                //
                // Maximum velocity constraint
                //

                if (midDisplacement[n][m] > maxDisplacement) {
                    if (verbosity > VerbosityNormal) {
                        cout << "Pairing old #" << oldCluster->id 
                                << " with new cluster #" << n 
                                << " rejected because of violation of max displacement."
                                << " (dR = " << midDisplacement[n][m]
                                << ", dR_max = " << maxDisplacement 
                                << endl;
                    }
                    continue;
                }

                typename Histogram<T>::ptr newHistogram;
                typename Histogram<T>::ptr oldHistogram;

//              #if WITH_OPENMP
//              #pragma omp critical
//              #endif
                {
                    newHistogram = newCluster->histogram(tracking_var_index,valid_min,valid_max);
                    oldHistogram = oldCluster->histogram(tracking_var_index,valid_min,valid_max);
 
                    // Processes in nature develop within certain bounds. It is 
                    // not possible that a cloud covers 10 pixels in one scan 
                    // and unit_multiplier in the next. The size deviation 
                    // constraint is created to prohibit matches between objects,
                    // which vary too much in size

                    // calculate cluster sizes by summing up all
                    // bins in their histograms.
                    size_t oldHistSum = (T) oldHistogram->sum();
                    size_t newHistSum = (T) newHistogram->sum();

                    // Determine the larger and smaller of the two
                    // numbers and their difference
                    T max_H = (T) max(oldHistSum, newHistSum);
                    T min_H = (T) min(oldHistSum, newHistSum);

                    histDiff[n][m] = (max_H == 0) ? 1.0 : (max_H - min_H);
                    
                    if (m_size_weight != 0.0) {
                        if ( histDiff[n][m] > maxHistD ) {
                            maxHistD = histDiff[n][m];
                        }
                    }
                    
                    // Now calculate the relative difference in percent
                    // (of number of points) between the two clusters: 
                    T size_deviation = (max_H - min_H) / min_H;
                        
                    // if the deviation is too big, exclude this match
                    if ( size_deviation > m_max_size_deviation ) {
                        cout << "Pairing old #" << oldCluster->id 
                                << " with new cluster #" << n 
                                << " rejected because of violation of max histogram size restraint:"
                                << " (dH = " << size_deviation
                                << ", dH_max = " << m_max_size_deviation
                                << endl;

                        continue;
                    }

                    // Just as the number of values will have some continuity,
                    // so will the distribution of values for a cluster. This
                    // fact is checked by creating a correlation between the
                    // histograms of the two clusters. Perfect match means
                    // a value of 1. No correlation at all means a value of 0. 
                    if ( m_corr_weight != 0.0 ) {
                        rank_correlation[n][m] = newHistogram->correlate_kendall( oldHistogram );
                    }
                }
                
                // only if all constraints are passed, the flag is set to true
                constraints_satisified[n][m] = true;

                if (constraints_satisified[n][m]) {
                    // Keep track of the largest distance in all possible
                    // matches for making values relative later.
                    if (midDisplacement[n][m] > maxMidD) {
                        maxMidD = midDisplacement[n][m];
                    }
                }
            } // done finishing correlation table

            if ( verbosity >= VerbosityNormal ) {
                cout << " done. (" << stop_timer() << "s)" << endl;
            }

            // Can't have zeros here
            if (maxHistD==0) {
                maxHistD = 1;
            }
        
    #pragma mark -
    #pragma mark Calculate probabilities

            if ( verbosity >= VerbosityNormal ) {
                cout << "Calculating matching table ... " << flush;;
                start_timer();
            }

//            #if WITH_OPENMP
//            #pragma omp parallel for schedule(dynamic)
//            #endif
            for (size_t idx = 0; idx < mapping.size(); idx++)
            {
                vector<int> index_pair = mapping.linear_to_grid(idx);
                
                int n = index_pair[0];
                int m = index_pair[1];

                typename Cluster<T>::ptr newCluster = current->clusters[n];
                typename Cluster<T>::ptr oldCluster = previous->clusters[m];

                // Only calculate values for pairs, that satisfy the
                // overlap constraint

                if (constraints_satisified[n][m])
                {
                    // The probability of the distance based matching estimate
                    // is the 
                    float prob_r = erfc( midDisplacement[n][m].get() / maxMidD.get() );

                    float prob_h = 0;
                    float prob_t = 0;

                    if (tracking_var_index >= 0)
                    {
                        // The probability for the histogram match is 
                        // the complementary error function of the relative
                        // histogram difference. The larger it is, the smaller
                        // the value will be
                        prob_h = erfc(histDiff[n][m]/maxHistD);
                        
                        // The probability of the 'signature' match is the
                        // raw output of the kendall's tau correlation of 
                        // the cluster histograms
                        prob_t = rank_correlation[n][m];
                    }

                    // The final matching probability is a 
                    // weighed sum of all three factors. 
                    sum_prob[n][m] 
                            = m_dist_weight * prob_r 
                            + m_size_weight * prob_h 
                            + m_corr_weight * prob_t;
                }
            }

            // print correlation table
            
            if ( verbosity >= VerbosityAll ) {
                
                for (int n=0; n < new_count; n++) {
                    typename Cluster<T>::ptr newCluster = current->clusters[n];
                    cout << endl << "Correlating new Cluster #" << n
                    << " (histogram size " << newCluster->histogram(tracking_var_index,valid_min,valid_max)->sum() << ")"
                    << " at " << newCluster->mode << endl;
                    for (int m=0; m < old_count; m++) {
                        if (constraints_satisified[n][m] || coverOldByNew[n][m] > 0.0 || coverNewByOld[n][m] > 0.0) {
                            typename Cluster<T>::ptr oldCluster = previous->clusters[m];
                            printf("\t<ID#%4lu>:\t(|H|=%5lu)\t\tdR=%4.1f\tdH=%5.4f\ttau=%7.4f\tsum=%6.4f\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                                    oldCluster->id,
                                    oldCluster->histogram(tracking_var_index,valid_min,valid_max)->sum(),
                                    midDisplacement[n][m].get(),
                                    histDiff[n][m],
                                    rank_correlation[n][m],
                                    sum_prob[n][m],
                                    coverOldByNew[n][m],
                                    coverNewByOld[n][m]);
                        }
                    }
                }
            }
            
            if ( verbosity >= VerbosityNormal ) {
                cout << " done. (" << stop_timer() << "s)" << endl;
            }
            
            // Prior to matchmaking, set all ids in the current set
            // back to NO_ID
            // index them, otherwise the area calculations will be nonsense
            for (size_t n=0; n < current->clusters.size(); n++) {
                current->clusters[n]->id = m3D::NO_ID;
            }

    #pragma mark -
    #pragma mark Matchmaking

            int n,m;
            ::units::values::meters_per_second velocitySum = ::units::values::meters_per_second(0);
            int velocityClusterCount = 0;
            float currentMaxProb = numeric_limits<float>::max();
            int maxIterations = new_count * old_count;
            int iterations = 0;
            set< typename Cluster<T>::ptr > used_clusters;

            if ( verbosity >= VerbosityNormal ) {
                cout << "Matchmaking ... " << flush;;
                start_timer();
                if ( verbosity >= VerbosityDetails ) {
                    cout << endl << "-- Matching Results --" << endl;
                }
            }

            // put the matches in a special data structure
            typedef pair<size_t,T>      match_t;
            typedef vector< match_t >   matchlist_t;
            matchlist_t matches;
            
            for (size_t idx = 0; idx < mapping.size(); idx++) {
                vector<int> index_pair = mapping.linear_to_grid(idx);
                int n = index_pair[0];
                int m = index_pair[1];
                if (!constraints_satisified[n][m]) continue;
                matches.push_back(match_t(idx,sum_prob[n][m]));
            }
            
            // sort the matches in descending order of probability
            sort(matches.begin(), 
                 matches.end(),  
                 boost::bind(&match_t::second, _1) > boost::bind(&match_t::second, _2));

            // now simply use this list to figure matches out
            
            // keep track of the new clusters as well (#355)
            id_set_t matched_ids;

            // in case we're re-running the tracking, clear some properties
            current->tracked_ids.clear();
            current->merges.clear();
            current->splits.clear();
            current->new_ids.clear();
            current->dropped_ids.clear();
            
            if ( verbosity >= VerbosityDetails ) {
                cout << "Examining prospective matches:" << endl;
            }
            
            for (size_t mi=0; mi < matches.size(); mi++) {
                match_t match = matches.at(mi);
                
                // back to n/m indexes
                vector<int> pairing = mapping.linear_to_grid(match.first);
                int n = pairing[0];
                int m = pairing[1];
                
                typename Cluster<T>::ptr new_cluster = current->clusters[n];
                typename Cluster<T>::ptr old_cluster = previous->clusters[m];
                
                if ( verbosity >= VerbosityDetails ) {
                    printf("\tprob=%6.4f \told cluster #%lu \tnew cluster #%d : ",
                            match.second, old_cluster->id, n);
                }
                
                // if the old cluster was matched already in a match
                // with higher probability, then skip this one
                id_set_t::iterator fi = current->tracked_ids.find(old_cluster->id);
                if ( fi != current->tracked_ids.end()) {
                    if ( verbosity >= VerbosityDetails ) {
                        cout << "rejected because previous cluster #" 
                                << old_cluster->id  << " was tracked already"
                                << endl;
                    }
                    continue;
                }                

                // Or if the NEW cluster was matched earlier, skip as well 
                fi = matched_ids.find(new_cluster->id);
                if (fi != matched_ids.end()) {
                    if ( verbosity >= VerbosityDetails ) {
                        cout << "rejected because new cluster #" 
                                << n  << " was tracked already"
                                << endl;
                    }
                    continue;
                }
                
                // Getting here means, that the match is accepted.
               
                // ID is continued
                new_cluster->id = old_cluster->id;
                
                // calculate displacement
                new_cluster->displacement = new_cluster->geometrical_center() 
                        - old_cluster->geometrical_center();
                
                // Remove the previous cluster from the race
                used_clusters.insert(old_cluster);
                
                // Remember that this ID was tracked
                current->tracked_ids.insert(old_cluster->id);

                // remove the new cluster from the race
                matched_ids.insert(new_cluster->id);

                // Update for mean velocity calculation
                ::units::values::meters_per_second velocity = midDisplacement[n][m] / this->m_deltaT;
                velocitySum += velocity;
                velocityClusterCount++;
                
                if (verbosity >= VerbosityNormal) {
                    cout << "accepted. (velocity " << velocity.get() << " m/s)"
                         << endl;
                }
            }
            
            if ( verbosity >= VerbosityNormal ) {
                cout << " done. (" << stop_timer() << "s)" << endl;
            }

    #pragma mark -
    #pragma mark Tag unmatched clusters

            if ( verbosity >= VerbosityDetails )
                printf("\n-- New ID assignments --\n");

            // for all new clusters that have not been assigned 
            // an id from one of the previous clusters, pick a 
            // fresh one now
            for ( n=0; n < new_count; n++ ) {
                typename Cluster<T>::ptr c = current->clusters[n];
                if ( c->id == m3D::NO_ID ) 
                {
                    // grab a new ID and memorize it for lookup
                    c->id = next_id(highest_id);
                    current->new_ids.insert( c->id );
                    
                    if ( verbosity >= VerbosityNormal ) {
                        cout << "#" << n << " at " << c->mode
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
            
            if ( verbosity >= VerbosityDetails ) {
                cout << "\n-- Goners --" << endl;
            }

            for ( ci = previous->clusters.begin(); ci != previous->clusters.end(); ci++ ) {
                typename Cluster<T>::ptr c = *ci;

                typename set< typename Cluster<T>::ptr >::iterator f = used_clusters.find( c );
                if (f != used_clusters.end()) continue;

                if ( verbosity >= VerbosityDetails ) {
                    cout << "#" << c->id << " at " << c->mode << " |H|=" 
                            << c->histogram(tracking_var_index,valid_min,valid_max)->sum() 
                            << endl;
                }

                current->dropped_ids.insert(c->id);
                hadGoners = true;
            }
            
            if (!hadGoners && verbosity >= VerbosityDetails) {
                cout << "none." << endl;
            }

    #pragma mark -
    #pragma mark Merging

            if ( verbosity >= VerbosityNormal ) {
                cout << "Handing merges and splits ..." << flush;;
                start_timer();
                if ( verbosity >= VerbosityDetails ) {  
                    cout << endl << "-- Merges --" << endl;
                }
            }

            // list of IDs that are continued through the tracking
            id_set_t continued_merged_ids;
            
            // list of IDs used previously in merges
            id_set_t merged_cluster_ids;
            
            bool had_merges = false;
            
            for ( n=0; n < current->clusters.size(); n++) 
            {
                vector<int> candidates;
                typename Cluster<T>::ptr new_cluster = current->clusters[n];
                
                // keep track of the largest candidate, that is: the cluster from
                // the previous series, that covers most of the area of the current
                // cluster

                size_t maxSize = 0;
                size_t largestCandidateIndex = 0;

                for ( m=0; m < previous->clusters.size(); m++ ) 
                {
                    T percentCovered = coverOldByNew[n][m];
                    
                    if (percentCovered >= this->m_ms_threshold ) 
                    {
                        typename Cluster<T>::ptr c = previous->clusters[m];

                        // check if this candidate was already used in a 
                        // merging earlier
                        id_set_t::iterator fi = merged_cluster_ids.find(c->id);
                        if (fi != merged_cluster_ids.end()) {
                            continue;
                        }
                        
//                        fi = tracked_ids.find(c->id);
//                        if (fi != tracked_ids.end() {
//                            continue;
//                        }

                        // add to the list of candidates 
                        candidates.push_back(m);

                        // track the largest candidate
                        if (c->size() > maxSize) {
                            maxSize = c->size();
                            largestCandidateIndex = m;
                        }
                    }
                }

                if ( candidates.size() > 1 )  {
                    
                    had_merges = true;

                    // collect merged cluster ids
                    id_set_t merged_from;

                    for ( int i=0; i < candidates.size(); i++ ) {
                        typename Cluster<T>::ptr c = previous->clusters[candidates[i]];
                        merged_cluster_ids.insert(c->id);
                        merged_from.insert(c->id);
                    }
                    if ( verbosity >= VerbosityNormal )
                        cout << "clusters " << merged_from << " seem to have merged into cluster " << new_cluster->id << endl;
                    
                    // store in attributes
                    current->merges[new_cluster->id] = merged_from;

                    // store the given ID
                    size_t the_id = new_cluster->id;

                    // check if the ID was new?
                    id_set_t::const_iterator it = current->new_ids.find(the_id);

                    if (it == current->new_ids.end())
                    {
                        if ( coverNewByOld[n][largestCandidateIndex] > this->m_msc_threshold && m_continueIDs)
                        {
                            // if the biggest congruence is at least 75%, continue the track
                            typename Cluster<T>::ptr c = previous->clusters[largestCandidateIndex];
                            new_cluster->id = c->id;
                            current->tracked_ids.insert(new_cluster->id);
                            continued_merged_ids.insert(new_cluster->id);

                            if ( verbosity >= VerbosityNormal )
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

                            if ( verbosity >= VerbosityNormal )
                                printf("\ttrack ID#%lu ends. new track ID#%lu begins\n", the_id, new_cluster->id );
                        }
                    }
                    else
                    {
                        if ( verbosity >= VerbosityNormal )
                            printf("\ttrack ID#%lu continues\n", new_cluster->id );
                    }
                }
            }

            if ( ! had_merges )
                if ( verbosity >= VerbosityDetails )
                    cout << "none." << endl;

    #pragma mark -
    #pragma mark Splitting

            if ( verbosity >= VerbosityDetails )
                cout << endl << "-- Splits --" << endl;

            bool had_splits = false;

            for ( m=0; m < previous->clusters.size(); m++ )
            {
                typename Cluster<T>::ptr old_cluster = previous->clusters[m];

                // check if the id was already continued in a merge
                id_set_t::iterator fi = continued_merged_ids.find(old_cluster->id);
                if (fi != continued_merged_ids.end()) continue;
                
                // collect candidates
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
                            // add to split candidates and record the candidate
                            // with the most coverage
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

                    if ( verbosity >= VerbosityNormal )
                        cout << "cluster " << old_cluster->id << " seems to have split into clusters " << split_into << endl;

                    for ( int i=0; i < candidates.size(); i++ )
                    {
                        typename Cluster<T>::ptr c = current->clusters[ candidates[i] ];

                        // check if the tracked id's contains c->id
                        m3D::id_t the_id = c->id;
                        id_set_t::iterator it = current->tracked_ids.find(the_id);

                        // If the largest one of the new clusters is at least 75% of the
                        // size of new previous cluster, the ID of the previous cluster
                        // is continued in the largest candidate
                        bool continueID = (candidates[i] == largestCandidateIndex) && (maxCover > this->m_msc_threshold);
                        if (continueID && m_continueIDs)
                        {
                            c->id = old_cluster->id;
                            current->tracked_ids.insert(c->id);
                            if ( verbosity >= VerbosityNormal )
                                printf("\ttrack ID#%lu continues\n", c->id );
                        }
                        else
                        {
                            if ( it != current->tracked_ids.end() )
                            {
                                // remove from tracked id's, re-tag and add to new IDs
                                current->tracked_ids.erase(it);
                                c->id = next_id(highest_id);
                                if ( verbosity >= VerbosityNormal )
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
                cout << " done. (" << stop_timer() << "s)" << endl;
            }
            
            std::cout.precision(prec);

        }

        current->tracking_performed = true;
        current->highest_id = highest_id;
    }
}
    
#endif
