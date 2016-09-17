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

#include <algorithm>
#include <boost/bind.hpp>
#include <math.h>
#include <netcdf>
#include <set>
#include <utility>
#include <vector>

#include "tracking.h"

namespace m3D {
    using namespace ::units;
    using namespace utils;
    using namespace utils::vectors;

#pragma mark -
#pragma mark Defaults

    template<typename T>
    tracking_param_t
    Tracking<T>::defaultParams() {
        tracking_param_t params;
        params.range_weight = 1.0;
        params.size_weight = 1.0;
        params.correlation_weight = 0.0;
        params.continueIDs = true;
        params.mergeSplitThreshold = 1.0f / 3.0f;
        params.mergeSplitContinuationThreshold = 0.75;
        params.maxVelocity = values::meters_per_second(50.0);
        params.max_deltaT = values::s(915);
        params.useOverlapConstraint = true;
        params.max_size_deviation = 3;
        params.useMeanVelocityConstraint = false;
        params.meanVelocityPercentage = 0.0;
        params.verbosity = VerbosityNormal;
        params.tracking_variable = "__default__";
        params.write_vtk = false;
        return params;
    }

#pragma mark -
#pragma mark Experimental

    template<typename T>
    void
    Tracking<T>::advectClusters(typename Tracking<T>::tracking_run_t &run) {
        using namespace utils::vectors;

        for (size_t i = 0; i < run.M; i++) {
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

#pragma mark -
#pragma mark Preliminaries

    template<typename T>
    bool
    Tracking<T>::initialise(typename Tracking<T>::tracking_run_t &run) {
        bool skip_tracking = false;
        bool logDetails = m_params.verbosity >= VerbosityDetails;
        run.cs = NULL;

        if (logDetails) {
            cout << endl;
        }
        // Check if there is anything to do
        if (run.N == 0) {
            if (logDetails) {
                cout << "\tNo current clusters. Skipping tracking." << endl;
            }
            return true;
        }

        // Check if the feature variables match
        if (run.previous->variables != run.current->variables) {
            cerr << "FATAL:Incompatible feature variables in the cluster files:"
                 << " previous:" << run.previous->variables
                 << " current:" << run.current->variables
                 << endl;
            exit(EXIT_FAILURE);
        }

        run.N = run.current->clusters.size();
        run.M = run.previous->clusters.size();

        // Find out the highest id from the previous list
        run.highestId = run.previous->highest_id;
        if (run.highestId == NO_ID || (run.highestId == 0 && run.M > 0)) {
            // figure it out from the highest found ID
            run.highestId = 0;
            for (size_t ci = 0; ci < run.M; ci++) {
                if (run.previous->clusters[ci]->id > run.highestId) {
                    run.highestId = run.previous->clusters[ci]->id;
                }
            }
        }

        // Find out the highest uuid from the previous list
        run.highestUuid = run.previous->highest_uuid;
        if (run.highestUuid == NO_UUID || (run.highestUuid == 0 && run.M > 0)) {
            // figure it out from the highest found ID
            run.highestUuid = 0;
            for (size_t ci = 0; ci < run.M; ci++) {
                if (run.previous->clusters[ci]->id > run.highestUuid) {
                    run.highestUuid = run.previous->clusters[ci]->id;
                }
            }
        }
        if (logDetails) {
            cout << "\thighest used id is " << run.highestId << endl;
            cout << "\thighest used uuid is " << run.highestUuid << endl;
        }

        // Get us a coordinate system
        NcFile *infoFile = run.current->file == NULL
                           ? run.previous->file : run.current->file;

        if (!skip_tracking) {

            run.cs = new CoordinateSystem<T>(infoFile,
                                             run.current->dimensions,
                                             run.current->dimension_variables);

            // Check cluster sizes
            for (size_t i = 0; i < run.M; i++) {
                typename Cluster<T>::ptr c = run.previous->clusters[i];
                if (c->size() == 0) {
                    cerr << "ERROR: previous cluster " << c->id << " has no points!" << endl;
                }
            }
            for (size_t i = 0; i < run.N; i++) {
                typename Cluster<T>::ptr c = run.current->clusters[i];
                if (c->size() == 0) {
                    cerr << "ERROR: current cluster " << c->id << " has no points!" << endl;
                }
            }

            // figure out tracking variable index (for histogram)
            run.tracking_var_index = -1;
            run.haveHistogramInfo = false;
            if (m_params.tracking_variable == "__default__") {
                std::string variable = run.current->variables[0];
                m_params.tracking_variable = variable;
                run.haveHistogramInfo = true;
                if (logDetails) {
                    cout << "Choosing default variable for histogram correlation: "
                         << m_params.tracking_variable << endl;
                }
            } else {
                bool found_tracking_var = false;
                for (size_t i = 0; i < run.current->variables.size(); i++) {
                    if (run.current->variables[i] == m_params.tracking_variable) {
                        found_tracking_var = true;
                        run.tracking_var_index = i;
                        run.haveHistogramInfo = true;
                    }
                }
                if (!found_tracking_var) {
                    cerr << "FATAL:tracking variable " << m_params.tracking_variable
                         << " is not part of the feature variables" << endl;
                    exit(EXIT_FAILURE);
                }
            }

            // Valid min/max of tracking variable
            utils::netcdf::unpacked_limits(infoFile->getVar(m_params.tracking_variable), run.valid_min, run.valid_max);

            // check time difference and determine displacement restraints
            ::units::values::s p_time = ::units::values::s(run.previous->timestamp);
            ::units::values::s c_time = ::units::values::s(run.current->timestamp);
            if (p_time >= c_time) {
                cerr << "ERROR:previous timestamp not greater than current timestamp!" << endl;
                return true;
            }
            run.deltaT = c_time - p_time;

            if (run.deltaT > m_params.max_deltaT) {
                cerr << "ERROR:files too far apart in time. Time difference:"
                     << run.deltaT << "s, longest accepted:"
                     << m_params.max_deltaT << "s" << endl;
                return true;
            }

            // calculate max displacement based on time difference
            // and max velocity
            run.current->tracking_time_difference = run.deltaT.get();
            run.maxDisplacement = m_params.maxVelocity * run.deltaT;

            if (logDetails) {
                cout << "\tmax velocity constraint:" << m_params.maxVelocity
                     << " and dT:" << run.deltaT
                     << " results in dR_max:" << run.maxDisplacement
                     << endl;
            }

            // TODO: mean velocity constraint

            // Minimum object radius for overlap constraint
            // TODO: this point should to be replaced with velocity vectors from SMV/RMV estimates
            run.overlap_constraint_velocity = m_params.maxVelocity;
            run.overlap_constraint_radius = 0.5 * run.deltaT * run.overlap_constraint_velocity;
            if (logDetails) {
                cout << "\toverlap constraint velocity:" << run.overlap_constraint_velocity
                     << " overlap constraint radius:" << run.overlap_constraint_radius
                     << endl;
            }

            // Provide a mapping for matchmaking and parallelisation
            vector<size_t> index_dims;
            index_dims.push_back(run.N);
            index_dims.push_back(run.M);
            run.mapping = LinearIndexMapping(index_dims);

            // Bestow the current cluster list with fresh uuids
            ClusterUtils<T>::provideUuids(run.current, run.highestUuid);

            if (m_params.verbosity > VerbosityNormal) {
                cout << endl << "-- previous clusters --" << endl;
                run.previous->print();
                cout << endl << "-- current clusters --" << endl;
                run.current->print();
            }
        }

        return skip_tracking;
    }

#pragma mark -
#pragma mark Matching

    template<typename T>
    void
    Tracking<T>::calculateCorrelationData(typename Tracking<T>::tracking_run_t &run) {
        bool logDetails = m_params.verbosity >= VerbosityAll;
        if (logDetails) {
            cout << endl;
        }

        // Allocate the correlation data
        run.rankCorrelation = SimpleMatrix<T>::create_matrix(run.N, run.M);
        run.midDisplacement = SimpleMatrix<::units::values::m>::create_matrix(run.N, run.M);
        run.sizeDifference = SimpleMatrix<T>::create_matrix(run.N, run.M);
        run.coverOldByNew = SimpleMatrix<T>::create_matrix(run.N, run.M);
        run.coverNewByOld = SimpleMatrix<T>::create_matrix(run.N, run.M);
        run.matchPossible = SimpleMatrix<T>::create_flag_matrix(run.N, run.M);
        run.likelihood = SimpleMatrix<T>::create_matrix(run.N, run.M);

        run.maxSizeDifference = numeric_limits<int>::min();
        run.maxMidDisplacement = ::units::values::m(numeric_limits<T>::min());

        // Index the cluster lists for quick overlap calculations

        // need to label the clusters with a unique id before we
        // index them, otherwise the area calculations will be nonsense
        for (size_t n = 0; n < run.N; n++) {
            run.current->clusters[n]->id = n;
        }
        ClusterIndex<T> index_c(run.current->clusters, run.cs->get_dimension_sizes());
        ClusterIndex<T> index_p(run.previous->clusters, run.cs->get_dimension_sizes());

        // Employ a mapping to parallelize

        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);

            int n = index_pair[0];
            int m = index_pair[1];

            // set all constraint flags to false to start with
            run.matchPossible[n][m] = false;

            typename Cluster<T>::ptr c = run.current->clusters[n];
            typename Cluster<T>::ptr p = run.previous->clusters[m];
            if (logDetails) {
                cout << "\tmatch uuid:" << p->uuid
                     << " id:" << p->id
                     << " with uuid:" << c->uuid << " ";
            }

            // How many percent of old cluster's points are shared by
            // the new cluster?
            run.coverOldByNew[n][m] = index_c.occupation_ratio(p, c);

            // How many percent of the new cluster's points are shared
            // by the old cluster
            run.coverNewByOld[n][m] = index_p.occupation_ratio(c, p);

            // Calculate this for merge/splits
            vector<T> oldCenter, newCenter;
            vector<T> dx = c->geometrical_center() - p->geometrical_center();
            run.midDisplacement[n][m] = ::units::values::m(vector_norm(run.cs->to_meters(dx)));

            //
            // Overlap constraint
            //
            // if the object is so big, that overlap is required at the given advection velocity
            // then check if that is the case. If no overlap exists, prohibit the match by setting
            // the constraint to false. If no overlap is required, the constraint is simply set to
            // true, thus allowing a match.

            if (m_params.useOverlapConstraint) {
                ::units::values::m radius = p->radius(run.cs);
                bool requires_overlap = (radius >= run.overlap_constraint_radius);
                if (requires_overlap && run.coverOldByNew[n][m] == 0.0) {
                    if (logDetails) {
                        cout << "precluded: violation of overlap constraint." << endl;
                    }
                    continue;
                }
            }

            //
            // Growth/shrink rate constraint
            //
            // Processes in nature develop within certain bounds. It is
            // not possible that a cloud covers 10 pixels in one scan
            // and 10.000 in the next. The size deviation constraint is
            // created to prohibit matches between objects, which vary
            // too much in size
            T maxSize = (T) max(p->size(), c->size());
            T minSize = (T) min(p->size(), c->size());
            T sizeDiff = (maxSize == 0) ? 1.0 : (maxSize - minSize);
            run.sizeDifference[n][m] = sizeDiff;
            if (m_params.size_weight != 0.0) {
                if (run.sizeDifference[n][m] > run.maxSizeDifference) {
                    run.maxSizeDifference = run.sizeDifference[n][m];
                }
            }
            T sizeDeviation = (maxSize - minSize) / minSize;
            if (sizeDeviation > m_params.max_size_deviation) {
                if (logDetails) {
                    cout << "precluded: violation of size constraint"
                         << " (dH:" << sizeDeviation << " values"
                         << " ,dH_max:" << m_params.max_size_deviation << " values)."
                         << endl;
                }
                continue;
            }

            //
            // Maximum velocity constraint
            //

            if (run.midDisplacement[n][m] > run.maxDisplacement) {
                if (logDetails) {
                    cout << "precluded: violation of max displacement"
                         << " (dR:" << run.midDisplacement[n][m]
                         << " dR_max:" << run.maxDisplacement
                         << ")." << endl;
                }
                continue;
            }

            //
            // Histogram correlation values
            //

            // Just as the number of values will have some continuity,
            // so will the distribution of values for a cluster. This
            // fact is checked by creating a correlation between the
            // histograms of the two clusters. Perfect match means
            // a value of 1. No correlation at all means a value of 0.

            if (run.haveHistogramInfo && m_params.correlation_weight != 0.0) {
                typename Histogram<T>::ptr hist_p
                        = p->histogram(run.tracking_var_index, run.valid_min, run.valid_max);
                typename Histogram<T>::ptr hist_c
                        = c->histogram(run.tracking_var_index, run.valid_min, run.valid_max);
                run.rankCorrelation[n][m] = hist_c->correlate_kendall(hist_p);
            }

            //
            // only if all constraints are passed, the flag is set to true
            //
            run.matchPossible[n][m] = true;

            // Keep track of the largest distance in all possible
            // matches for making values relative later.
            if (run.midDisplacement[n][m] > run.maxMidDisplacement) {
                run.maxMidDisplacement = run.midDisplacement[n][m];
            }

            if (logDetails) {
                cout << "possible." << endl;
            }

        } // done finishing correlation table

        // Now re-tag new clusters with no id
        run.current->erase_identifiers();

        // Can't have zeros here
        if (run.maxSizeDifference == 0) {
            run.maxSizeDifference = 1;
        }
    }

    template<typename T>
    void
    Tracking<T>::calculateProbabilities(typename Tracking<T>::tracking_run_t &run) {
        bool logDetails = m_params.verbosity >= VerbosityAll;
        if (logDetails) {
            cout << endl;
        }

        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);

            int n = index_pair[0];
            int m = index_pair[1];

            typename Cluster<T>::ptr c = run.current->clusters[n];
            typename Cluster<T>::ptr p = run.previous->clusters[m];

            // Only calculate values for pairs, that satisfy the
            // overlap constraint

            if (run.matchPossible[n][m]) {
                // The probability of the distance based matching estimate
                // is the 
                float prob_r = erfc(run.midDisplacement[n][m].get() / run.maxMidDisplacement.get());

                float prob_h = 0;
                float prob_t = 0;

                if (run.tracking_var_index >= 0) {
                    // The probability for the histogram match is 
                    // the complementary error function of the relative
                    // histogram difference. The larger it is, the smaller
                    // the value will be
                    prob_h = erfc(run.sizeDifference[n][m] / run.maxSizeDifference);

                    // The probability of the 'signature' match is the
                    // raw output of the kendall's tau correlation of 
                    // the cluster histograms
                    prob_t = run.rankCorrelation[n][m];
                }

                // The final matching probability is a 
                // weighed sum of all three factors. 
                run.likelihood[n][m]
                        = m_params.range_weight * prob_r
                          + m_params.size_weight * prob_h
                          + m_params.correlation_weight * prob_t;
            }
        }

        // print correlation table

        if (m_params.verbosity >= VerbosityAll) {

            for (int n = 0; n < run.N; n++) {
                typename Cluster<T>::ptr c = run.current->clusters[n];

                if (logDetails) {
                    size_t histSize = c->size();
                    cout << "\tcluster uuid:" << c->uuid << " ("
                         << " |H|:" << histSize << ")"
                         << " mode:" << c->mode
                         << ")" << endl;
                }

                for (int m = 0; m < run.M; m++) {
                    if (run.matchPossible[n][m]) {
                        typename Cluster<T>::ptr p = run.previous->clusters[m];
                        printf("\t\tuuid:%4llu \tid:%4lu\t(|H|=%5lu)\t\tdR=%4.1f\tdH=%5.4f\ttau=%7.4f\tsum=%6.4f\t\tcovON=%3.2f\t\tcovNO=%3.2f\n",
                               p->uuid,
                               p->id,
                               p->size(),
                               run.midDisplacement[n][m].get(),
                               run.sizeDifference[n][m],
                               run.rankCorrelation[n][m],
                               run.likelihood[n][m],
                               run.coverOldByNew[n][m],
                               run.coverNewByOld[n][m]);
                    }
                }
            }
        }
    }

    template<typename T>
    void
    Tracking<T>::matchmaking(typename Tracking<T>::tracking_run_t &run) {
        bool logAll = m_params.verbosity >= VerbosityAll;
        bool logDetails = m_params.verbosity >= VerbosityDetails;
        if (logDetails) {
            cout << endl;
        }

        int n, m;
        ::units::values::meters_per_second velocitySum = ::units::values::meters_per_second(0);
        int velocityClusterCount = 0;
        float currentMaxProb = numeric_limits<float>::max();
        int maxIterations = run.N * run.M;
        int iterations = 0;

        // put the matches in a special data structure
        for (size_t idx = 0; idx < run.mapping.size(); idx++) {
            vector<int> index_pair = run.mapping.linear_to_grid(idx);
            int n = index_pair[0];
            int m = index_pair[1];
            if (!run.matchPossible[n][m]) continue;
            run.matches.push_back(match_t(idx, run.likelihood[n][m]));
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

        if (logDetails) {
            cout << "\tmatches:" << endl;
        }

        for (size_t mi = 0; mi < run.matches.size(); mi++) {
            match_t match = run.matches.at(mi);

            // back to n/m indexes
            vector<int> pairing = run.mapping.linear_to_grid(match.first);
            int n = pairing[0];
            int m = pairing[1];

            typename Cluster<T>::ptr c = run.current->clusters[n];
            typename Cluster<T>::ptr p = run.previous->clusters[m];

            // if the old cluster was matched already in a match
            // with higher probability, then skip this one
            id_set_t::iterator fi = run.current->tracked_ids.find(p->id);
            if (fi != run.current->tracked_ids.end()) {
//                if (logAll) {
//                    cout << "\t\tuuid:" << old_cluster->uuid << " id:" << old_cluster->id
//                         << " with uuid:" << new_cluster->uuid << " ";
//                    cout << "rejected: previous cluster already matched." << endl;
//                }
                continue;
            }

            // Or if the NEW cluster was matched earlier, skip as well 
            fi = run.matched_uuids.find(c->uuid);
            if (fi != run.matched_uuids.end()) {
//                if (logAll) {
//                    cout << "\t\tuuid:" << old_cluster->uuid << " id:" << old_cluster->id
//                         << " with uuid:" << new_cluster->uuid << " ";
//                    cout << "rejected: new cluster already matched." << endl;
//                }
                continue;
            }

            // Getting here means, that the match is accepted.

            // ID is continued
            c->id = p->id;

            // calculate displacement
            c->displacement = c->geometrical_center()
                              - p->geometrical_center();

            // Remember that this ID was tracked, removing
            // the old cluster from the race
            run.current->tracked_ids.insert(p->id);

            // remove the new cluster from the race
            run.matched_uuids.insert(c->uuid);

            // Update for mean velocity calculation
            ::units::values::meters_per_second velocity = run.midDisplacement[n][m] / run.deltaT;
            velocitySum += velocity;
            velocityClusterCount++;

            if (logDetails) {
                cout << "\t\tuuid:" << p->uuid << " id:" << p->id
                     << " with uuid:" << c->uuid << " ("
                     << " accepted: (velocity: " << velocity
                     << " displacement: " << c->displacement << ")." << endl;
            }
        }

        run.averageVelocity = velocitySum / boost::numeric_cast<T>(velocityClusterCount);

        // Tag unmatched clusters
        bool hadFreshIds = false;
        if (logDetails) {
            cout << "\tnew track ids:" << endl;
        }
        // for all new clusters that have not been assigned
        // an id from one of the previous clusters, pick a 
        // fresh one now
        for (int n = 0; n < run.N; n++) {
            typename Cluster<T>::ptr c = run.current->clusters[n];
            if (c->id == m3D::NO_ID) {
                // Get a fresh ID
                c->id = nextId(run.highestId);
                if (logDetails) {
                    cout << "\t\tuuid: " << c->uuid << " id:" << c->id << "." << endl;
                }
                run.current->new_ids.insert(c->id);
                hadFreshIds = true;
            }
        }
        if (logDetails && !hadFreshIds) {
            cout << "\t\tnone." << endl;
        }

        //Determine drop-outs
        bool hadGoners = false;
        typename Cluster<T>::list::iterator ci;
        if (logDetails) {
            cout << "\tterminated track ids:" << endl;
        }
        for (ci = run.previous->clusters.begin(); ci != run.previous->clusters.end(); ci++) {
            typename Cluster<T>::ptr c = *ci;
            id_set_t::const_iterator f = run.current->tracked_ids.find(c->id);
            if (f != run.current->tracked_ids.end()) continue;
            if (logDetails) {
                cout << "\t\tid: " << c->id << endl;
            }
            run.current->dropped_ids.insert(c->id);
            hadGoners = true;
        }
        if (!hadGoners && logDetails) {
            cout << "\t\tnone." << endl;
        }
    }

#pragma mark -
#pragma mark Merging

    template<typename T>
    double
    Tracking<T>::getMergeCriterion(typename Tracking<T>::tracking_run_t &run, const int &n, const int &m) {
        typename Cluster<T>::ptr c = run.current->clusters[n];
        typename Cluster<T>::ptr p = run.previous->clusters[m];
        double nbo = run.coverNewByOld[n][m];
        double dR = run.midDisplacement[n][m].get();
        double dH = run.sizeDifference[n][m];
        double s = nbo * (erfc(dR / run.maxMidDisplacement.get()) + erfc(dH / run.maxSizeDifference));
        return s;
    }

    template<typename T>
    int
    Tracking<T>::findBestMergeCandidate(typename Tracking<T>::tracking_run_t &run,
                                        const int &n,
                                        const vector<int> &candidates) {
        double maxS = -1.0;
        int maxM = -1;
        bool maxIsTied = false;
        for (size_t i = 0; i < candidates.size(); i++) {
            int m = candidates[i];
            if (run.coverNewByOld[n][m] >= m_params.mergeSplitContinuationThreshold) {
                double s = getMergeCriterion(run, n, m);
                if (s >= maxS) {
                    if (s == maxS) maxIsTied = true;
                    maxM = m;
                    maxS = s;
                }
            }
        }
        return maxIsTied ? -1 : maxM;
    }

    template<typename T>
    void
    Tracking<T>::getMergeCandidates(typename Tracking<T>::tracking_run_t &run,
                                    const int &n,
                                    bool &track_flag,
                                    vector<int> &candidates,
                                    id_set_t &candidateIds) {
        typename Cluster<T>::ptr c = run.current->clusters.at(n);
        for (size_t m = 0; m < run.M; m++) {
            typename Cluster<T>::ptr p = run.previous->clusters.at(m);
            if (c->id == p->id) {
                candidates.push_back(m);
                candidateIds.insert(p->id);
                track_flag = true;
            } else {

                // Check if this (previous) cluster was already associated
                // with a cluster from the current set? If so, exclude it
                id_set_t::const_iterator fi = run.current->tracked_ids.find(p->id);
                if (fi != run.current->tracked_ids.end()) continue;

                // Not tracked
                T obn = run.coverOldByNew[n][m];
                if (obn >= m_params.mergeSplitThreshold) {

                    // Check for each of the candidates what the merge criteria
                    // s is and find out if there is another combination with
                    // higher s. If so, do not add the candidate here.
                    double s1 = getMergeCriterion(run, n, m);
                    bool foundBetter = false;
                    for (int nn = 0; nn < run.N && !foundBetter; nn++) {
                        if (nn == n) continue;
                        double s2 = getMergeCriterion(run, nn, m);
                        if (s2 > s1) {
                            foundBetter = true;
                        }
                    }
                    if (!foundBetter) {
                        candidateIds.insert(p->id);
                        candidates.push_back(m);
                    }
                }
            }
        }
    }

    template<typename T>
    void
    Tracking<T>::handleMerges(typename Tracking<T>::tracking_run_t &run) {
        bool logDetails = m_params.verbosity >= VerbosityDetails;
        if (logDetails) {
            cout << endl << "\t-- Merges" << endl;
        }
        bool had_merges = false;
        size_t n, m;

        for (n = 0; n < run.N; n++) {
            typename Cluster<T>::ptr c = run.current->clusters[n];

            bool track_flag = false;
            id_set_t candidateIds;
            vector<int> candidates;
            getMergeCandidates(run, n, track_flag, candidates, candidateIds);

            if (candidates.size() > 1) {
                had_merges = true;
                if (logDetails) {
                    cout << "\t\tprevious ids:" << candidateIds << " have merged into "
                         << "uuid:" << c->uuid << " id:" << c->id << "." << endl;
                }
                if (!track_flag && m_params.continueIDs) {
                    int winner = findBestMergeCandidate(run, n, candidates);
                    if (winner >= 0) {
                        typename Cluster<T>::ptr p = run.previous->clusters.at(winner);
                        if (logDetails) {
                            cout << "\t\t\tuuid:" << c->uuid << " id:" << c->id
                                 << " re-tagged with id:" << p->id
                                 << endl
                                 << "\t\t\ttrack id:" << p->id << " continues." << endl;
                        }
                        run.current->new_ids.erase(c->id);
                        run.current->tracked_ids.insert(p->id);
                        c->id = p->id;
                    }
                }
                id_set_t merged_ids;
                for (int i = 0; i < candidates.size(); i++) {
                    typename Cluster<T>::ptr p = run.previous->clusters.at(candidates[i]);
                    merged_ids.insert(p->id);
                }
                run.current->merges[c->id] = merged_ids;
            }
        }

        if (logDetails && !had_merges) {
            cout << "\t\tnone." << endl;
        }
    }

#pragma mark -
#pragma mark Splitting

    template<typename T>
    double
    Tracking<T>::getSplitCriterion(typename Tracking<T>::tracking_run_t &run, const int &n, const int &m) {
        typename Cluster<T>::ptr c = run.current->clusters[n];
        typename Cluster<T>::ptr p = run.previous->clusters[m];
        double obn = run.coverOldByNew[n][m];
        double dR = run.midDisplacement[n][m].get();
        double dH = run.sizeDifference[n][m];
        double s = obn * (erfc(dR / run.maxMidDisplacement.get()) + erfc(dH / run.maxSizeDifference));
        return s;
    }

    template<typename T>
    int
    Tracking<T>::findBestSplitCandidate(typename Tracking<T>::tracking_run_t &run,
                                        const int &m,
                                        const vector<int> &candidates) {
        double maxS = -1.0;
        int maxM = -1;
        bool maxIsTied = false;
        for (size_t i = 0; i < candidates.size(); i++) {
            int n = candidates[i];
            if (run.coverOldByNew[n][m] >= m_params.mergeSplitContinuationThreshold) {
                double s = getSplitCriterion(run, n, m);
                if (s >= maxS) {
                    if (s == maxS) maxIsTied = true;
                    maxM = n;
                    maxS = s;
                }
            }
        }
        return maxIsTied ? -1 : maxM;
    }

    template<typename T>
    void
    Tracking<T>::getSplitCandidates(typename Tracking<T>::tracking_run_t &run,
                                    const int &m,
                                    bool &track_flag,
                                    vector<int> &candidates,
                                    uuid_set_t &candidateUuids) {
        typename Cluster<T>::ptr p = run.previous->clusters.at(m);
        for (size_t n = 0; n < run.N; n++) {
            typename Cluster<T>::ptr c = run.current->clusters.at(n);
            if (c->id == p->id) {
                candidates.push_back(n);
                candidateUuids.insert(c->uuid);
                track_flag = true;
            } else {

                // Check if this (current) cluster was already associated
                // with a cluster from the previous set? If so, exclude it
                id_set_t::const_iterator fi = run.current->tracked_ids.find(c->id);
                if (fi != run.current->tracked_ids.end()) continue;

                // Not tracked
                T obn = run.coverNewByOld[n][m];
                if (obn >= m_params.mergeSplitThreshold) {

                    // Check for each of the candidates what the split criteria
                    // s is and find out if there is another combination with
                    // higher s. If so, do not add the candidate here.

                    double s1 = getSplitCriterion(run, n, m);
                    bool foundBetter = false;
                    for (int mm = 0; mm < run.M && !foundBetter; mm++) {
                        if (mm == m) continue;
                        double s2 = getSplitCriterion(run, n, mm);
                        if (s2 > s1) {
                            foundBetter = true;
                        }
                    }
                    if (!foundBetter) {
                        candidateUuids.insert(c->uuid);
                        candidates.push_back(n);
                    }
                }
            }
        }
    }


    template<typename T>
    void
    Tracking<T>::handleSplits(typename Tracking<T>::tracking_run_t &run) {
        bool logDetails = m_params.verbosity >= VerbosityDetails;
        if (logDetails) {
            cout << "\t-- Splits" << endl;
        }
        bool had_splits = false;
        size_t n, m;
        for (m = 0; m < run.M; m++) {
            typename Cluster<T>::ptr p = run.previous->clusters[m];

            bool track_flag = false;
            uuid_set_t candidateUuids;
            vector<int> candidates;
            getSplitCandidates(run, m, track_flag, candidates, candidateUuids);

            if (candidates.size() > 1) {
                had_splits = true;
                if (logDetails) {
                    cout << "\t\tuuid:" << p->uuid << " id:" << p->id
                         << " split into uuids:" << candidateUuids << "." << endl;
                }
                if (!track_flag && m_params.continueIDs) {
                    int winner = findBestSplitCandidate(run, m, candidates);
                    if (winner >= 0) {
                        typename Cluster<T>::ptr c = run.current->clusters.at(winner);
                        if (logDetails) {
                            cout << "\t\t\tuuid:" << c->uuid << " id:" << c->id
                                 << " re-tagged with id:" << p->id
                                 << endl
                                 << "\t\t\ttrack id:" << p->id << " continues." << endl;
                        }
                        run.current->new_ids.erase(c->id);
                        run.current->tracked_ids.insert(c->id);
                        c->id = p->id;
                    }
                }

                // Compose the split_ids set
                id_set_t split_ids;
                for (int i = 0; i < candidates.size(); i++) {
                    typename Cluster<T>::ptr c = run.current->clusters.at(candidates[i]);
                    split_ids.insert(c->id);
                }
                run.current->splits[p->id] = split_ids;
            }
        }

        if (logDetails && !had_splits) {
            cout << "\t\tnone." << endl;
        }
    }

    template<typename T>
    void
    Tracking<T>::removeScheduled(const tracking_run_t &run) {
        id_set_t::const_iterator si;
        for (si = run.scheduled_for_removal.begin(); si != run.scheduled_for_removal.end(); ++si) {
            m3D::id_t id = *si;
            run.current->tracked_ids.erase(id);
        }
    }

#pragma mark -
#pragma mark Entry point

    template<typename T>
    void
    Tracking<T>::track(typename ClusterList<T>::ptr previous,
                       typename ClusterList<T>::ptr current) {
        bool logNormal = m_params.verbosity >= VerbosityNormal;

        if (logNormal) cout << "Tracking:" << endl;

        tracking_run_t run;
        run.previous = previous;
        run.current = current;
        if (logNormal) start_timer("-- Calculating preliminaries ... ");
        bool skip_tracking = initialise(run);
        if (logNormal) stop_timer("done");
        if (skip_tracking) {
            return;
        }

        // set precision for outputting floating point numbers
        std::streamsize prec = std::cout.precision();
        std::cout.precision(2);
        cout << setiosflags(ios::fixed);

        // Advect the previous clusters by their displacement
        // vectors to make the results more accurate
        // utils::VisitUtils<T>::write_clusters_vtu_wholesale(previous, cs, "unshifted");
        if (m_params.useDisplacementVectors) {
            if (logNormal) start_timer("-- Shifting clusters ... ");
            advectClusters(run);
            if (logNormal) stop_timer("done");

        }
        // utils::VisitUtils<T>::write_clusters_vtu_wholesale(previous, cs, "shifted");

        // Set all ids in the current set back to NO_ID
        current->erase_identifiers();

        if (logNormal) start_timer("-- Calculating correlation data ... ");
        calculateCorrelationData(run);
        if (logNormal) stop_timer("done");

        if (logNormal) start_timer("-- Calculating match probabilities ... ");
        calculateProbabilities(run);
        if (logNormal) stop_timer("done");

        if (logNormal) start_timer("-- Matchmaking ... ");
        matchmaking(run);
        if (logNormal) stop_timer("done");

        if (logNormal) start_timer("-- Merging and splitting ... ");
        // Calculate merges and splits
        handleMerges(run);
        handleSplits(run);
        removeScheduled(run);
        if (logNormal) stop_timer("done");

        // Restore save precision
        std::cout.precision(prec);

        current->tracking_performed = true;
        current->highest_id = run.highestId;
        current->highest_uuid = run.highestUuid;

        // Clean up
        if (run.cs != NULL) {
            delete run.cs;
            run.cs = NULL;
        }
    }
}


#endif
