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

#include <meanie3D/utils/vector_utils.h>
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

    using namespace utils;
    using namespace utils::vectors;

    template <typename T>
    void
    Tracking<T>::advectClusters(typename Tracking<T>::tracking_run_t &run)
    {
        using namespace utils::vectors;

        for (size_t i = 0; i < run.previous->clusters.size(); i++) {
            typename Cluster<T>::ptr c = run.previous->clusters[i];
            if (vector_norm<T>(c->displacement) != 0) {
                // shift all points by the displacement vector
                for (size_t i = 0; i < c->size(); i++) {
                    typename Point<T>::ptr p = c->get_points()[i];
                    // shift coordinate
                    p->coordinate += c->displacement;
                    // TODO: shift gridded coordinate
                    for (size_t ci = 0; ci < c->displacement.size(); ci++) {
                        // shift feature space point spatial range
                        p->values[ci] += c->displacement[ci];
                    }
                }
            }
        }
    }

    template <typename T>
    bool
    Tracking<T>::calculatePreliminaries(typename Tracking<T>::tracking_run_t &run)
    {

        bool skip_tracking = false;

        // Check if there is anything to do
        if (run.current->clusters.size() == 0) {
            if (m_verbosity >= VerbosityNormal) {
                cout << endl << "No current clusters. Skipping tracking." << endl;
            }
            return true;
        }

        // Find out the highest id from the previous list

        run.highestId = run.previous->highest_id;

        if (run.highestId == NO_ID || (run.highestId == 0 && run.previous->clusters.size() > 0)) {
            // figure it out from the highest found ID
            run.highestId = 0;
            for (size_t ci = 0; ci < run.previous->clusters.size(); ci++) {
                if (run.previous->clusters[ci]->id > run.highestId) {
                    run.highestId = run.previous->clusters[ci]->id;
                }
            }
        }

        if (m_verbosity >= VerbosityDetails) {
            cout << "Highest used ID is " << run.highestId << endl;
        }

        if (!skip_tracking) {

            // Check cluster sizes
            for (size_t i = 0; i < run.previous->clusters.size(); i++) {
                typename Cluster<T>::ptr c = run.previous->clusters[i];
                if (c->size() == 0) {
                    cerr << "ERROR: previous cluster " << c->id << " has no points!" << endl;
                }
            }
            for (size_t i = 0; i < run.current->clusters.size(); i++) {
                typename Cluster<T>::ptr c = run.current->clusters[i];
                if (c->size() == 0) {
                    cerr << "ERROR: current cluster " << c->id << " has no points!" << endl;
                }
            }

            // figure out tracking variable index (for histogram)
            run.tracking_var_index = -1;
            if (!run.tracking_variable.empty()) {
                for (size_t i = 0; i < run.current->feature_variables.size(); i++) {
                    if (run.current->feature_variables[i].getName() == run.tracking_variable) {
                        run.tracking_var_index = i;
                        break;
                    }
                }

                // TODO: fail only if histogram weight is > 0
                if (run.tracking_var_index < 0) {
                    cerr << "FATAL:Illegal variable name " << run.tracking_variable << " for histogram comparison " << endl;
                    exit(EXIT_FAILURE);
                }
                // Valid min/max of tracking variable
                utils::netcdf::get_valid_range(run.current->ncFile->getVar(run.tracking_variable), run.valid_min, run.valid_max);
            }

            // check time difference and determine displacement restraints
            ::units::values::s p_time = ::units::values::s(run.previous->timestamp);
            ::units::values::s c_time = ::units::values::s(run.current->timestamp);
            if (p_time >= c_time) {
                cerr << "ERROR:previous timestamp not greater than current timestamp!" << endl;
                return true;
            }
            run.deltaT = c_time - p_time;

            if (run.deltaT > m_max_deltaT) {
                cerr << "ERROR:files too far apart in time. Time difference:"
                        << run.deltaT << "s, longest accepted:"
                        << m_max_deltaT << "s" << endl;
                return true;
            }

            // calculate max displacement based on time difference
            // and max velocity
            run.current->tracking_time_difference = run.deltaT.get();
            run.maxDisplacement = this->m_maxVelocity * run.deltaT;

            if (m_verbosity >= VerbosityDetails) {
                cout << "Max velocity constraint at "
                        << this->m_maxVelocity.get() << "m/s"
                        << " and dT=" << run.deltaT.get()
                        << " results in dRmax=" << run.maxDisplacement.get()
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

            run.overlap_constraint_velocity = m_maxVelocity;
            run.overlap_constraint_radius = 0.5 * run.deltaT * run.overlap_constraint_velocity;

            // Provide a mapping
            vector<size_t> index_dims;
            index_dims.push_back(run.current->clusters.size());
            index_dims.push_back(run.previous->clusters.size());
            run.mapping = LinearIndexMapping(index_dims);
        }

        return skip_tracking;
    }

    template <typename T>
    void
    Tracking<T>::calculateCorrelationData(typename Tracking<T>::tracking_run_t &run)
    {

        size_t old_count = run.previous->clusters.size();
        size_t new_count = run.current->clusters.size();

        // Allocate the correlation data
        run.rank_correlation = SimpleMatrix<T>::create_matrix(new_count, old_count);
        run.midDisplacement = SimpleMatrix< ::units::values::m >::create_matrix(new_count, old_count);
        run.histDiff = SimpleMatrix<T>::create_matrix(new_count, old_count);
        run.sum_prob = SimpleMatrix<T>::create_matrix(new_count, old_count);
        run.coverOldByNew = SimpleMatrix<T>::create_matrix(new_count, old_count);
        run.coverNewByOld = SimpleMatrix<T>::create_matrix(new_count, old_count);
        run.constraints_satisified = SimpleMatrix<T>::create_flag_matrix(new_count, old_count);

        run.maxHistD = numeric_limits<int>::min();
        run.maxMidD = ::units::values::m(numeric_limits<T>::min());

        // Index the cluster lists for quick overlap calculations

        // need to label the clusters with a unique id before we
        // index them, otherwise the area calculations will be nonsense
        for (size_t n = 0; n < run.current->clusters.size(); n++) {
            run.current->clusters[n]->id = n;
        }
        ClusterIndex<T> new_index(run.current->clusters, run.cs->get_dimension_sizes());
        ClusterIndex<T> old_index(run.previous->clusters, run.cs->get_dimension_sizes());

        // Employ a mapping to parallelize


        //            #if WITH_OPENMP
        //            #pragma omp parallel for schedule(dynamic)
        //            #endif
        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);

            int n = index_pair[0];
            int m = index_pair[1];

            // set all constraint flags to false to start with
            run.constraints_satisified[n][m] = false;

            typename Cluster<T>::ptr newCluster = run.current->clusters[n];
            typename Cluster<T>::ptr oldCluster = run.previous->clusters[m];

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
                T old_covered_by_new = new_index.occupation_ratio(oldCluster, newCluster);
                run.coverOldByNew[n][m] = old_covered_by_new;

                // How many percent of the new cluster's points are shared
                // by the old cluster
                T new_covered_by_old = old_index.occupation_ratio(newCluster, oldCluster);
                run.coverNewByOld[n][m] = new_covered_by_old;

                if (m_useOverlapConstraint) {
                    ::units::values::m radius = oldCluster->radius(run.cs);
                    bool requires_overlap = (radius >= run.overlap_constraint_radius);
                    if (requires_overlap) {
                        overlap_constraint_satisfied = run.coverOldByNew[n][m] > 0.0;
                    }
                }

                if (overlap_constraint_satisfied) {
                    // calculate average mid displacement, that is the 
                    // distance between the two clusters geometrical
                    // centers. 

                    vector<T> oldCenter, newCenter;
                    oldCenter = oldCluster->geometrical_center();
                    newCenter = newCluster->geometrical_center();

                    vector<T> dx = newCenter - oldCenter;
                    run.midDisplacement[n][m] = ::units::values::m(vector_norm(run.cs->to_meters(dx)));
                }
            }

            if (!overlap_constraint_satisfied && m_verbosity > VerbosityNormal) {
                cout << "Pairing old #" << oldCluster->id
                        << " with new cluster #" << n
                        << " rejected because of violation of overlap constraint"
                        << endl;
                continue;
            }

            //
            // Maximum velocity constraint
            //

            if (run.midDisplacement[n][m] > run.maxDisplacement) {
                if (m_verbosity > VerbosityNormal) {
                    cout << "Pairing old #" << oldCluster->id
                            << " with new cluster #" << n
                            << " rejected because of violation of max displacement."
                            << " (dR = " << run.midDisplacement[n][m]
                            << ", dR_max = " << run.maxDisplacement
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
                newHistogram = newCluster->histogram(run.tracking_var_index, run.valid_min, run.valid_max);
                oldHistogram = oldCluster->histogram(run.tracking_var_index, run.valid_min, run.valid_max);

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

                run.histDiff[n][m] = (max_H == 0) ? 1.0 : (max_H - min_H);

                if (m_size_weight != 0.0) {
                    if (run.histDiff[n][m] > run.maxHistD) {
                        run.maxHistD = run.histDiff[n][m];
                    }
                }

                // Now calculate the relative difference in percent
                // (of number of points) between the two clusters: 
                T size_deviation = (max_H - min_H) / min_H;

                // if the deviation is too big, exclude this match
                if (size_deviation > m_max_size_deviation) {
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
                if (m_corr_weight != 0.0) {
                    run.rank_correlation[n][m] = newHistogram->correlate_kendall(oldHistogram);
                }
            }

            // only if all constraints are passed, the flag is set to true
            run.constraints_satisified[n][m] = true;

            if (run.constraints_satisified[n][m]) {
                // Keep track of the largest distance in all possible
                // matches for making values relative later.
                if (run.midDisplacement[n][m] > run.maxMidD) {
                    run.maxMidD = run.midDisplacement[n][m];
                }
            }
        } // done finishing correlation table

        if (m_verbosity >= VerbosityNormal) {
            cout << " done. (" << stop_timer() << "s)" << endl;
        }

        // Can't have zeros here
        if (run.maxHistD == 0) {
            run.maxHistD = 1;
        }
    }

#pragma mark -
#pragma mark Calculate probabilities

    template <typename T>
    void
    Tracking<T>::calculateProbabilities(typename Tracking<T>::tracking_run_t & run)
    {
        size_t old_count = run.previous->clusters.size();
        size_t new_count = run.current->clusters.size();

        //            #if WITH_OPENMP
        //            #pragma omp parallel for schedule(dynamic)
        //            #endif
        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);

            int n = index_pair[0];
            int m = index_pair[1];

            typename Cluster<T>::ptr newCluster = run.current->clusters[n];
            typename Cluster<T>::ptr oldCluster = run.previous->clusters[m];

            // Only calculate values for pairs, that satisfy the
            // overlap constraint

            if (run.constraints_satisified[n][m]) {
                // The probability of the distance based matching estimate
                // is the 
                float prob_r = erfc(run.midDisplacement[n][m].get() / run.maxMidD.get());

                float prob_h = 0;
                float prob_t = 0;

                if (run.tracking_var_index >= 0) {
                    // The probability for the histogram match is 
                    // the complementary error function of the relative
                    // histogram difference. The larger it is, the smaller
                    // the value will be
                    prob_h = erfc(run.histDiff[n][m] / run.maxHistD);

                    // The probability of the 'signature' match is the
                    // raw output of the kendall's tau correlation of 
                    // the cluster histograms
                    prob_t = run.rank_correlation[n][m];
                }

                // The final matching probability is a 
                // weighed sum of all three factors. 
                run.sum_prob[n][m]
                        = m_dist_weight * prob_r
                        + m_size_weight * prob_h
                        + m_corr_weight * prob_t;
            }
        }

        // print correlation table

        if (m_verbosity >= VerbosityAll) {

            for (int n = 0; n < new_count; n++) {
                typename Cluster<T>::ptr newCluster = run.current->clusters[n];
                cout << endl << "Correlating new Cluster #" << n
                        << " (histogram size " << newCluster->histogram(run.tracking_var_index, run.valid_min, run.valid_max)->sum() << ")"
                        << " at " << newCluster->mode << endl;
                for (int m = 0; m < old_count; m++) {
                    if (run.constraints_satisified[n][m] || run.coverOldByNew[n][m] > 0.0 || run.coverNewByOld[n][m] > 0.0) {
                        typename Cluster<T>::ptr oldCluster = run.previous->clusters[m];
                        printf("\t<ID#%4lu>:\t(|H|=%5lu)\t\tdR=%4.1f\tdH=%5.4f\ttau=%7.4f\tsum=%6.4f\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                                oldCluster->id,
                                oldCluster->histogram(run.tracking_var_index, run.valid_min, run.valid_max)->sum(),
                                run.midDisplacement[n][m].get(),
                                run.histDiff[n][m],
                                run.rank_correlation[n][m],
                                run.sum_prob[n][m],
                                run.coverOldByNew[n][m],
                                run.coverNewByOld[n][m]);
                    }
                }
            }
        }

        // Prior to matchmaking, set all ids in the current set
        // back to NO_ID
        // index them, otherwise the area calculations will be nonsense
        for (size_t n = 0; n < run.current->clusters.size(); n++) {
            run.current->clusters[n]->id = m3D::NO_ID;
        }
    }

    template <typename T>
    void
    Tracking<T>::matchmaking(typename Tracking<T>::tracking_run_t & run)
    {
        size_t old_count = run.previous->clusters.size();
        size_t new_count = run.current->clusters.size();

        int n, m;
        ::units::values::meters_per_second velocitySum = ::units::values::meters_per_second(0);
        int velocityClusterCount = 0;
        float currentMaxProb = numeric_limits<float>::max();
        int maxIterations = new_count * old_count;
        int iterations = 0;

        if (m_verbosity >= VerbosityDetails) {
            cout << endl << "-- Matching Results --" << endl;
        }

        // put the matches in a special data structure
        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);
            int n = index_pair[0];
            int m = index_pair[1];
            if (!run.constraints_satisified[n][m]) continue;
            run.matches.push_back(match_t(idx, run.sum_prob[n][m]));
        }

        // sort the matches in descending order of probability
        sort(run.matches.begin(),
                run.matches.end(),
                boost::bind(&match_t::second, _1) > boost::bind(&match_t::second, _2));

        // in case we're re-running the tracking, clear some properties
        run.current->tracked_ids.clear();
        run.current->merges.clear();
        run.current->splits.clear();
        run.current->new_ids.clear();
        run.current->dropped_ids.clear();

        if (m_verbosity >= VerbosityDetails) {
            cout << "Examining prospective matches:" << endl;
        }

        for (size_t mi = 0; mi < run.matches.size(); mi++) {
            match_t match = run.matches.at(mi);

            // back to n/m indexes
            vector<int> pairing = run.mapping.linear_to_grid(match.first);
            int n = pairing[0];
            int m = pairing[1];

            typename Cluster<T>::ptr new_cluster = run.current->clusters[n];
            typename Cluster<T>::ptr old_cluster = run.previous->clusters[m];

            if (m_verbosity >= VerbosityNormal) {
                cout << "Pairing old #" << old_cluster->id
                        << " with new cluster #" << new_cluster->id
                        << " ";
            }

            // if the old cluster was matched already in a match
            // with higher probability, then skip this one
            id_set_t::iterator fi = run.current->tracked_ids.find(old_cluster->id);
            if (fi != run.current->tracked_ids.end()) {
                if (m_verbosity >= VerbosityNormal) {
                    cout << "rejected because previous cluster #"
                            << old_cluster->id << " was tracked already"
                            << endl;
                }
                continue;
            }

            // Or if the NEW cluster was matched earlier, skip as well 
            fi = run.matched_ids.find(new_cluster->id);
            if (fi != run.matched_ids.end()) {
                if (m_verbosity >= VerbosityNormal) {
                    cout << "rejected because new cluster #"
                            << n << " was tracked already"
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
            run.used_clusters.insert(old_cluster);

            // Remember that this ID was tracked
            run.current->tracked_ids.insert(old_cluster->id);

            // remove the new cluster from the race
            run.matched_ids.insert(new_cluster->id);

            // Update for mean velocity calculation
            ::units::values::meters_per_second velocity = run.midDisplacement[n][m] / run.deltaT;
            velocitySum += velocity;
            velocityClusterCount++;

            if (m_verbosity >= VerbosityNormal) {
                cout << "accepted. (velocity " << velocity.get() << " m/s"
                        << ",displacement=" << new_cluster->displacement << ")"
                        << endl;
            }
        }

        run.averageVelocity = velocitySum / boost::numeric_cast<T>(velocityClusterCount);

        // Tag unmatched clusters

        if (m_verbosity >= VerbosityDetails)
            printf("\n-- New ID assignments --\n");

        // for all new clusters that have not been assigned 
        // an id from one of the previous clusters, pick a 
        // fresh one now
        for (int n = 0; n < run.current->clusters.size(); n++) {
            typename Cluster<T>::ptr c = run.current->clusters[n];
            if (c->id == m3D::NO_ID) {
                // grab a new ID and memorize it for lookup
                c->id = next_id(run.highestId);
                run.current->new_ids.insert(c->id);

                if (m_verbosity >= VerbosityNormal) {
                    cout << "#" << n << " at " << c->mode
                            << " with |H|=" << c->histogram(run.tracking_var_index, run.valid_min, run.valid_max)->sum()
                            << " is assigned new ID #" << c->id
                            << endl;
                }
            }
        }

        //Determine drop-outs

        bool hadGoners = false;
        typename Cluster<T>::list::iterator ci;

        if (m_verbosity >= VerbosityDetails) {
            cout << "\n-- Goners --" << endl;
        }

        for (ci = run.previous->clusters.begin(); ci != run.previous->clusters.end(); ci++) {
            typename Cluster<T>::ptr c = *ci;

            typename set< typename Cluster<T>::ptr >::iterator f = run.used_clusters.find(c);
            if (f != run.used_clusters.end()) continue;

            if (m_verbosity >= VerbosityDetails) {
                cout << "#" << c->id << " at " << c->mode << " |H|="
                        << c->histogram(run.tracking_var_index, run.valid_min, run.valid_max)->sum()
                        << endl;
            }

            run.current->dropped_ids.insert(c->id);
            hadGoners = true;
        }

        if (!hadGoners && m_verbosity >= VerbosityDetails) {
            cout << "none." << endl;
        }
    }

    template <typename T>
    void
    Tracking<T>::handleMerges(typename Tracking<T>::tracking_run_t &run)
    {
        if (m_verbosity >= VerbosityDetails) {
            cout << endl << "-- Merges --" << endl;
        }

        bool had_merges = false;

        size_t n, m;

        for (n = 0; n < run.current->clusters.size(); n++) {
            vector<int> candidates;
            typename Cluster<T>::ptr new_cluster = run.current->clusters[n];

            // keep track of the largest candidate, that is: the cluster from
            // the previous series, that covers most of the area of the current
            // cluster

            size_t maxSize = 0;
            size_t largestCandidateIndex = 0;

            for (m = 0; m < run.previous->clusters.size(); m++) {
                T percentCovered = run.coverOldByNew[n][m];

                if (percentCovered >= this->m_ms_threshold) {
                    typename Cluster<T>::ptr c = run.previous->clusters[m];

                    // check if this candidate was already used in a 
                    // merging earlier
                    id_set_t::iterator fi = run.merged_cluster_ids.find(c->id);
                    if (fi != run.merged_cluster_ids.end()) {
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

            if (candidates.size() > 1) {

                had_merges = true;

                // collect merged cluster ids
                id_set_t merged_from;

                for (int i = 0; i < candidates.size(); i++) {
                    typename Cluster<T>::ptr c = run.previous->clusters[candidates[i]];
                    run.merged_cluster_ids.insert(c->id);
                    merged_from.insert(c->id);
                }
                if (m_verbosity >= VerbosityNormal)
                    cout << "clusters " << merged_from << " seem to have merged into cluster " << new_cluster->id << endl;

                // store in attributes
                run.current->merges[new_cluster->id] = merged_from;

                // store the given ID
                size_t the_id = new_cluster->id;

                // check if the ID was new?
                id_set_t::const_iterator it = run.current->new_ids.find(the_id);

                if (it == run.current->new_ids.end()) {
                    if (run.coverNewByOld[n][largestCandidateIndex] > this->m_msc_threshold && m_continueIDs) {
                        // if the biggest congruence is at least 75%, continue the track
                        typename Cluster<T>::ptr c = run.previous->clusters[largestCandidateIndex];
                        new_cluster->id = c->id;
                        run.current->tracked_ids.insert(new_cluster->id);
                        run.continued_merged_ids.insert(new_cluster->id);

                        if (m_verbosity >= VerbosityNormal)
                            printf("\ttrack ID#%lu continues\n", new_cluster->id);
                    } else {
                        // Otherwise dish out a fresh tag and add to new IDs
                        new_cluster->id = next_id(run.highestId);
                        run.current->new_ids.insert(new_cluster->id);

                        // as well as remove from tracked IDs
                        id_set_t::iterator fi = run.current->tracked_ids.find(the_id);
                        if (fi != run.current->tracked_ids.end()) {
                            run.current->tracked_ids.erase(the_id);
                        }

                        if (m_verbosity >= VerbosityNormal)
                            printf("\ttrack ID#%lu ends. new track ID#%lu begins\n", the_id, new_cluster->id);
                    }
                } else {
                    if (m_verbosity >= VerbosityNormal)
                        printf("\ttrack ID#%lu continues\n", new_cluster->id);
                }
            }
        }

        if (m_verbosity >= VerbosityDetails && !had_merges)
            cout << "none." << endl;
    }

    template <typename T>
    void
    Tracking<T>::handleSplits(typename Tracking<T>::tracking_run_t &run)
    {
        if (m_verbosity >= VerbosityDetails)
            cout << endl << "-- Splits --" << endl;

        bool had_splits = false;
        size_t n, m;
        for (m = 0; m < run.previous->clusters.size(); m++) {
            typename Cluster<T>::ptr old_cluster = run.previous->clusters[m];

            // check if the id was already continued in a merge
            id_set_t::iterator fi = run.continued_merged_ids.find(old_cluster->id);
            if (fi != run.continued_merged_ids.end()) continue;

            // collect candidates
            vector<float> candidates;
            T maxCover = 0.0;
            size_t largestCandidateIndex = 0;

            for (n = 0; n < run.current->clusters.size(); n++) {
                if (run.constraints_satisified[n][m]) {
                    T overlap = run.coverNewByOld[n][m];
                    if (overlap > this->m_ms_threshold) {
                        // add to split candidates and record the candidate
                        // with the most coverage
                        candidates.push_back(n);

                        if (run.coverOldByNew[n][m] > maxCover) {
                            maxCover = run.coverOldByNew[n][m];
                            largestCandidateIndex = n;
                        }
                    }
                }
            }

            if (candidates.size() > 1) {
                had_splits = true;
                id_set_t split_into;

                for (int i = 0; i < candidates.size(); i++)
                    split_into.insert(run.current->clusters[candidates[i]]->id);

                if (m_verbosity >= VerbosityNormal)
                    cout << "cluster " << old_cluster->id << " seems to have split into clusters " << split_into << endl;

                for (int i = 0; i < candidates.size(); i++) {
                    typename Cluster<T>::ptr c = run.current->clusters[ candidates[i] ];

                    // check if the tracked id's contains c->id
                    m3D::id_t the_id = c->id;
                    id_set_t::iterator it = run.current->tracked_ids.find(the_id);

                    // If the largest one of the new clusters is at least 75% of the
                    // size of new previous cluster, the ID of the previous cluster
                    // is continued in the largest candidate
                    bool continueID = (candidates[i] == largestCandidateIndex) && (maxCover > this->m_msc_threshold);
                    if (continueID && m_continueIDs) {
                        c->id = old_cluster->id;
                        run.current->tracked_ids.insert(c->id);
                        if (m_verbosity >= VerbosityNormal)
                            printf("\ttrack ID#%lu continues\n", c->id);
                    } else {
                        if (it != run.current->tracked_ids.end()) {
                            // remove from tracked id's, re-tag and add to new IDs
                            run.current->tracked_ids.erase(it);
                            c->id = next_id(run.highestId);
                            if (m_verbosity >= VerbosityNormal)
                                printf("\tnew track ID#%lu begins\n", c->id);
                            run.current->new_ids.insert(c->id);
                        }
                    }
                }

                run.current->splits[old_cluster->id] = split_into;
            }
        }

        if (m_verbosity >= VerbosityDetails && !had_splits)
            cout << "none." << endl;
    }

    template <typename T>
    void
    Tracking<T>::track(typename ClusterList<T>::ptr previous,
            typename ClusterList<T>::ptr current,
            const CoordinateSystem<T> *cs,
            const std::string *tracking_variable_name)
    {
#pragma mark -
#pragma mark Preliminaries        
        if (m_verbosity >= VerbosityNormal) cout << "Tracking:" << endl;

        tracking_run_t run;
        run.previous = previous;
        run.current = current;
        run.cs = cs;

        bool skip_tracking = calculatePreliminaries(run);
        if (skip_tracking) {
            return;
        }

        // set precision for outputting correlation tables
        std::streamsize prec = std::cout.precision();
        std::cout.precision(2);
        cout << setiosflags(ios::fixed);

        // Advect the previous clusters by their displacement
        // vectors to make the results more accurate
        // utils::VisitUtils<T>::write_clusters_vtu_wholesale(previous, cs, "unshifted");
        if (this->m_useDisplacementVectors) {
            if (m_verbosity >= VerbosityNormal) start_timer("-- Shifting clusters ... ");
            this->advectClusters(run);
            if (m_verbosity >= VerbosityNormal) stop_timer("done");

        }
        // utils::VisitUtils<T>::write_clusters_vtu_wholesale(previous, cs, "shifted");

        if (m_verbosity >= VerbosityNormal) start_timer("-- Calculating correlation data ... ");
        calculateCorrelationData(run);
        if (m_verbosity >= VerbosityNormal) stop_timer("done");

        if (m_verbosity >= VerbosityNormal) start_timer("-- Calculating match probabilities ... ");
        calculateProbabilities(run);
        if (m_verbosity >= VerbosityNormal) stop_timer("done");

        if (m_verbosity >= VerbosityNormal) start_timer("-- Matchmaking ... ");
        matchmaking(run);
        if (m_verbosity >= VerbosityNormal) stop_timer("done");

        if (m_verbosity >= VerbosityNormal) start_timer("-- Merging and splitting ... ");
        handleMerges(run);
        handleSplits(run);
        if (m_verbosity >= VerbosityNormal) stop_timer("done");

        // Restore save precision
        std::cout.precision(prec);

        current->tracking_performed = true;
        current->highest_id = run.highestId;
    }
}


#endif
