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
#include <meanie3D/utils/matrix.h>
#include <meanie3D/utils.h>
#include <meanie3D/clustering.h>

#include <netcdf>

#include "tracking.h"

namespace m3D {

    using namespace utils;

    /** This class contains the tracking code.
     */
    template <typename T>
    class Tracking
    {
    private:

        // Member Variables

        T m_dist_weight; // correlation weight distance
        T m_size_weight; // correlation weight histogram sum
        T m_corr_weight; // correlation weight histogram rank correlation

        ::units::values::s m_max_deltaT; // What is the maximum time between slices for valid tracking (in seconds)
        float m_max_size_deviation; // how many percent may the objects vary in size (number of points) between scans?

        ::units::values::meters_per_second m_maxVelocity; // physical maximum speed of objects in m/s
        bool m_useOverlapConstraint; // make use of max velocity constraint?

        bool m_useMeanVelocityConstraint; // make use of max velocity constraint?
        T m_meanVelocitySecurityPercentage; // Percentual amount of deviation from mean velocity allowed

        float m_ms_threshold; // Percentage of coverage required to be a candidate for merge/split?
        float m_msc_threshold; // How many percent of coverage is required in splits/merges for track continuation?

        bool m_continueIDs; // Continue ids through merging/splitting?
        bool m_useDisplacementVectors; // Experimental: use displacement vectors to advect previous clusters?

        Verbosity m_verbosity; // Verbosity of output

    public:

        /** Constructor
         * @param weight for distance correlation
         * @param weight for size correlation
         * @param weight for histogram rank correlation
         */
        Tracking(T wr = 1.0, T wd = 1.0, T wt = 1.0,
                ::units::values::s max_delta_t = ::units::values::s(930),
                const Verbosity verbosity = VerbosityNormal)
        : m_dist_weight(wr)
        , m_size_weight(wd)
        , m_corr_weight(wt)
        , m_max_deltaT(max_delta_t) // 15 minutes )(plus 30 seconds slack)
        , m_useMeanVelocityConstraint(false) // limit deviation from mean velocity (false)
        , m_meanVelocitySecurityPercentage(0.5) // to 50 %
        , m_maxVelocity(::units::values::meters_per_second(100.0)) // limit max velocity to 30 m/s (~108 km/h)
        , m_max_size_deviation(2.5) // how many percent may the objects vary in size (number of points) between scans (250%)
        , m_useOverlapConstraint(true)
        , m_ms_threshold(0.5) // percentage coverage old/new for merge/split (33%)
        , m_msc_threshold(0.75) // percentage coverage old/new in merge/split for continuing track (75%)
        , m_continueIDs(true) // continue IDs through splits/merges ?)
        , m_useDisplacementVectors(false)
        , m_verbosity(verbosity)
        {
        };

        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param clusters form the current run
         * @param coordinate system
         * @param name of variable to be used for the histogram correlation. 
         *        If NULL, histogram comparison is ignored completely.
         * @param verbosity
         */
        void track(typename ClusterList<T>::ptr previous,
                typename ClusterList<T>::ptr current,
                const CoordinateSystem<T> *cs,
                const std::string *tracking_variable_name = NULL);

        // Accessors

        /** Maximum allowed speed for object pairings. 
         * 
         * @return 
         */
        ::units::values::meters_per_second maxTrackingSpeed()
        {
            return m_maxVelocity;
        }

        void setMaxTrackingSpeed(::units::values::meters_per_second speed)
        {
            m_maxVelocity = speed;
        }

        /** How much area is covered between clusters from previous and
         * current time slice to make it a split/merge? 
         * 
         * @return 
         */
        float mergeSplitThreshold()
        {
            return m_ms_threshold;
        }

        void setMergeSplitThreshold(float value)
        {
            m_ms_threshold = value;
        }

        /** What is the minimum coverage between clusters from the previous
         * and current time slice after a merge / split occured to consider
         * the continuation of the ID?
         * 
         * @return 
         */
        float mergeSplitContinuationThreshold()
        {
            return m_msc_threshold;
        }

        void setMergeSplitContinuationThreshold(float value)
        {
            m_msc_threshold = value;
        }

        /** Are IDs even considered to be continued in merges/splits?
         * 
         * @return 
         */
        bool continueIDs()
        {
            return m_continueIDs;
        }

        void setContinueIDs(bool value)
        {
            m_continueIDs = value;
        }

        /** 
         * Are displacement vectors used to improve tracking?
         * @return 
         */
        bool useDisplacementVectors()
        {
            return m_useDisplacementVectors;
        }

        void setUseDisplacementVectors(bool value)
        {
            m_useDisplacementVectors = value;
        }

        /**
         * 
         * @return 
         */
        ::units::values::s maxDeltaT()
        {
            return m_max_deltaT;
        }

        void setMaxDeltaT(::units::values::s seconds)
        {
            m_max_deltaT = seconds;
        }

        /**
         * 
         * @param verbosity
         */
        void setVerbosity(Verbosity verbosity)
        {
            m_verbosity = verbosity;
        }

    protected:

        typedef pair<size_t, T> match_t;
        typedef vector< match_t > matchlist_t;

        /**
         * Bundles data that constitutes a tracking run. These are mostly
         * things derived at the beginning, such as bounds, derived parameters
         * or similar. Also the correlation data between clusters.
         */
        typedef struct
        {
            typename ClusterList<T>::ptr current; // previous cluster list
            typename ClusterList<T>::ptr previous; // current cluster list

            // Properties
            const CoordinateSystem<T> *cs; // Coordinate system (for transformations)
            LinearIndexMapping mapping; // maps i <-> (n,m)

            m3D::id_t highestId;        // Stores the highest used ID
            m3D::uuid_t highestUuid;    // Stores the highest used UUID

            std::string tracking_variable;
            int tracking_var_index;
            T valid_min;
            T valid_max;

            ::units::values::s deltaT; // Difference in seconds
            ::units::values::m maxDisplacement; // Calculated from deltaT and max velocity
            ::units::values::m maxMidD;
            int maxHistD;

            ::units::values::m overlap_constraint_radius;
            ::units::values::meters_per_second overlap_constraint_velocity;

            id_set_t matched_uuids; // uuid of new clusters that were matched

            matchlist_t matches;    // final matching result

            id_set_t scheduled_for_removal; // Set of ids to be removed at the end of the run.

            ::units::values::meters_per_second averageVelocity;
            // Correlation data
            typename SimpleMatrix<T>::matrix_t rank_correlation;
            typename SimpleMatrix< ::units::values::m>::matrix_t midDisplacement;
            typename SimpleMatrix<T>::matrix_t histDiff;
            typename SimpleMatrix<T>::matrix_t sum_prob;
            typename SimpleMatrix<T>::matrix_t coverOldByNew;
            typename SimpleMatrix<T>::matrix_t coverNewByOld;
            typename SimpleMatrix<T>::flag_matrix_t constraints_satisified;

        } tracking_run_t;

        /**
         * Move the clusters in the list by their displacement vectors 
         * (which usually have been obtained on the last tracking run)
         * Only clusters with existing displacement vectors are moved. 
         * 
         * @param clusters
         */
        void advectClusters(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         * @return 
         */
        bool calculatePreliminaries(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void calculateCorrelationData(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void calculateProbabilities(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void matchmaking(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void handleMerges(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void handleSplits(typename Tracking<T>::tracking_run_t &run);

        /**
         * 
         * @param run
         */
        void removeScheduled(const tracking_run_t &run);

        /**
         *
         */
        bool tagIfNeeded(typename Tracking<T>::tracking_run_t &run,
                         typename Cluster<T>::ptr c);
    };
}

#endif
