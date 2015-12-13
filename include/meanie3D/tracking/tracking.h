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

#ifndef M3D_TRACKING_CLASS_H
#define M3D_TRACKING_CLASS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/array/linear_index_mapping.h>
#include <meanie3D/clustering/cluster.h>
#include <meanie3D/clustering/cluster_list.h>
#include <meanie3D/utils/matrix.h>
#include <meanie3D/utils/time_utils.h>

#include <netcdf>

#include "tracking.h"

namespace m3D {

    using namespace utils;

    /**
     * Parameters to the tracking algorithm.
     */
    typedef struct {
        std::string previous_filename; // Path to previous cluster file
        std::string current_filename; // Path to current cluster file

        bool write_vtk; // Write out clusters in .vtk? (TODO: remove)
        std::vector<size_t> vtk_dimension_indexes; // See --vtk-dimensions in meanie3D-detect

        double range_weight;        // Correlation weight for distance between clusters
        double size_weight;         // Correlation weight for size comparison
        double correlation_weight;  // Correlation weight for histogram comparison
        std::string tracking_variable; // Variable name to be used for histogram correlation.

        ::units::values::meters_per_second maxVelocity; // what is the allowed top speed of objects?
        ::units::values::s max_deltaT; // How much time in seconds allowed between scans?
        double max_size_deviation; // How many percent may the objects vary in size (number of points) between scans?

        bool continueIDs; // Continue ids across merges and splits?
        double mergeSplitThreshold; // What is the minimum overlap for merges/splits?
        double mergeSplitContinuationThreshold; // What is the overlap required for continuing ids in merges/splits?

        bool useDisplacementVectors; // Use displacment vectors (experimental)

        bool useOverlapConstraint; // Do objects require overlap to be tracked (if their speed/size is low/big enough)
        bool useMeanVelocityConstraint; // Are object matches constraint by average velocity? (currently defunct)
        double meanVelocityPercentage; // Maximum allowed deviation from average velocity? (currently defunct)

        Verbosity verbosity; // Output level (0 (none) to 3 (all))

    } tracking_param_t;


    /** This class contains the tracking code.
     */
    template <typename T>
    class Tracking {

    protected:

        typedef pair<size_t, T> match_t;
        typedef vector<match_t> matchlist_t;

        /**
         * Bundles data that constitutes a tracking run. These are mostly
         * things derived at the beginning, such as bounds, derived parameters
         * or similar. Also the correlation data between clusters.
         */
        typedef struct {
            // Properties
            typename ClusterList<T>::ptr current; // previous cluster list
            typename ClusterList<T>::ptr previous; // current cluster list
            size_t N,M;                     // Shortcuts for lenghts of previous and current lists.
            const CoordinateSystem<T> *cs;  // Coordinate system (for transformations)
            LinearIndexMapping mapping;     // maps i <-> (n,m)

            m3D::id_t highestId;        // Stores the highest used ID
            m3D::uuid_t highestUuid;    // Stores the highest used UUID

            bool haveHistogramInfo;             // Indicates if histogram based data is available
            int tracking_var_index;             // Index of the variable used for histogram
            T valid_min;                        // lower valid bound of histogram variable
            T valid_max;                        // upper valid bound of histogram variable

            ::units::values::s deltaT;          // Difference in seconds
            ::units::values::m maxDisplacement; // Calculated from deltaT and max velocity
            ::units::values::m maxMidDisplacement;         // Maximum center displacment
            int maxSizeDifference;                       // Maximum allowed histogram size difference

            ::units::values::meters_per_second averageVelocity;
            ::units::values::m overlap_constraint_radius;
            ::units::values::meters_per_second overlap_constraint_velocity;

            id_set_t matched_uuids;             // uuid of new clusters that were matched
            matchlist_t matches;                // final matching result
            id_set_t scheduled_for_removal;     // Set of ids to be removed at the end of the run.

            // Correlation data
            typename SimpleMatrix<T>::matrix_t rankCorrelation;
            typename SimpleMatrix< ::units::values::m >::matrix_t midDisplacement;
            typename SimpleMatrix<T>::matrix_t sizeDifference;
            typename SimpleMatrix<T>::matrix_t likelihood;
            typename SimpleMatrix<T>::matrix_t coverOldByNew;
            typename SimpleMatrix<T>::matrix_t coverNewByOld;
            typename SimpleMatrix<T>::flag_matrix_t matchPossible;

        } tracking_run_t;

    private:

        tracking_param_t m_params;

    public:

        /**
         * Obtain a set of default parameters for the tracking.
         * @return default params.
         */
        static tracking_param_t defaultParams();

        /**
         * Instantiates an instance with default params.
         */
        Tracking() : m_params(Tracking::defaultParams()) {};

        /** Constructor
         * @param weight for distance correlation
         * @param weight for size correlation
         * @param weight for histogram rank correlation
         */
        Tracking(tracking_param_t params) : m_params(params) {};

        /**
         * Runs the meanie3D tracking algorithm.
         * @param current
         * @param previous
         */
        void track(typename ClusterList<T>::ptr previous,
                   typename ClusterList<T>::ptr current);

    protected:

        /**
         * Experimental!!
         * Move the clusters in the list by their displacement vectors 
         * (which usually have been obtained on the last tracking run)
         * Only clusters with existing displacement vectors are moved. 
         * 
         * @param clusters
         */
        void advectClusters(typename Tracking<T>::tracking_run_t &run);

        /**
         * Initialises a tracking run.
         * @param run
         * @return 
         */
        bool initialise(typename Tracking<T>::tracking_run_t &run);

        /**
         * Calculates data to base match probabilities on.
         * @param run
         */
        void calculateCorrelationData(typename Tracking<T>::tracking_run_t &run);

        /**
         * Calculates match probabilities.
         * @param run
         */
        void calculateProbabilities(typename Tracking<T>::tracking_run_t &run);

        /**
         * Analyses the probabilities to find the most likely matches.
         * @param run
         */
        void matchmaking(typename Tracking<T>::tracking_run_t &run);

        /**
         * Called after matchmaking. Analyses and tags splits.
         * @param run
         */
        void handleSplits(typename Tracking<T>::tracking_run_t &run);

        /**
         * Called after matchmaking. Analyses and tags merges.
         * @param run
         */
        void handleMerges(typename Tracking<T>::tracking_run_t &run);

        /**
         *
         * @param run
         */
        void removeScheduled(const tracking_run_t &run);

        /**
         * Gets the split criterion value. The value is calculated considering:
         * <ul>
         *  <li>dR - geometrical center difference</li>
         *  <li>dH - size difference</li>
         *  <li>obn - overlap percentage (old by new)</li>
         * </ul>
         * Calculates: s = erfc(dR) + erfc(dH) + erf(obn)
         * @param tracking context
         * @param index of current cluster
         * @param index of previous cluster
         * @returns split criterion value
         */
        double
                getSplitCriterion(typename Tracking<T>::tracking_run_t &run,
                                  const int &m,
                                  const int &n);


        /**
         * Examine split criterion in unclear situations to find the
         * optimal candidate. Examines all possible candidates
         * by calculating the split criterion s:
         * =&gt; max(s) is the winner
         * =&gt; if equal candidates, no one wins
         * @param run tracking context
         * @param index of the cluster that is the result of the merge.
         * @param indexes of candidates
         * @return -1 if no one wins, number of candidate else.
         */
        int findBestSplitCandidate(typename Tracking<T>::tracking_run_t &run,
                                   const int &m,
                                   const vector<int> &candidates);


        /**
         *
         */
        void getSplitCandidates(typename Tracking<T>::tracking_run_t &run,
                                const int& m,
                                bool &track_flag,
                                vector<int> &candidates,
                                uuid_set_t &candidateUuids);

        /**
         * Examine merge criterion in unclear situations to find the
         * optimal candidate. Examines all possible candidates
         * by calculating the merge criterion s:
         * =&gt; max(s) is the winner
         * =&gt; if equal candidates, no one wins
         * @param run tracking context
         * @param index of the cluster that is the result of the merge.
         * @param indexes of candidates
         * @return -1 if no one wins, number of candidate else.
         */
        int findBestMergeCandidate(typename Tracking<T>::tracking_run_t &run,
                                   const int &n,
                                   const vector<int> &indexes);


        /**
         * Gets the merge criterion value. The value is calculated considering:
         * <ul>
         *  <li>dR - geometrical center difference</li>
         *  <li>dH - size difference</li>
         *  <li>nbo - overlap percentage (new by old)</li>
         * </ul>
         * Calculates: s = erfc(dR) + erfc(dH) + erf(nbo)
         * @param tracking context
         * @param index of current cluster
         * @param index of previous cluster
         * @returns merge criterion value.
         */
        double getMergeCriterion(typename Tracking<T>::tracking_run_t &run,
                                 const int &n,
                                 const int &m);

        /**
         * Compile a list of candidates that might have merged into c.
         * If there is a match, add it. Only consider other candidates
         * that have not been matched. Only consider those, if they
         * satisfy the overlap criterion for merge and split.
         *
         * @param run tracking context
         * @param index of merged cluster
         * @param track_flag contains <code>true</code> if the id of the merged
         * cluster was found in any of the candidates after the call.
         * @param contains indexes of candidates after the call.
         * @param contains ids of candidates after the call.
         */
        void getMergeCandidates(typename Tracking<T>::tracking_run_t &run,
                                const int& n,
                                bool &track_flag,
                                vector<int> &candidateIndexes,
                                id_set_t &candidateIds);

    };
}

#endif
