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
#ifndef M3D_DETECTION_H
#define    M3D_DETECTION_H

#include <map>
#include <vector>
#include <netcdf>
#include <string>

#include <meanie3D/filters/replacement_filter.h>
#include <meanie3D/filters/scalespace_filter.h>
#include <meanie3D/index/search_parameters.h>

namespace m3D {

    using namespace utils;
    using namespace std;
    using namespace netCDF;

    /**
     * Data structure for specifying a detection run. 
     */
    template<typename T>
    struct detection_params_t
    {
        // Input file
        string filename;

        // Name of resulting cluster file
        string output_filename;

        // Dimensions in the netCDF file to use.
        vector<std::string> dimensions;

        // Spatial range of Featurespace
        vector<std::string> dimension_variables;

        // Value range of featurespace 
        vector<std::string> variables;

        // When running detection in netCDF files with time dimension,
        // this needs to be set to the index in the time dimension to use.
        int time_index;

        // Maps index in feature-space variables to lower threshold. 
        map<int, double> lower_thresholds;

        // Maps index in feature-space variables to upper threshold. 
        map<int, double> upper_thresholds;

        // Map from index in feature-space variables to a value to replace
        // values outside of valid range or threshold with.
        map<int, double> replacement_values;

        // Minimum cluster size in number of grid points. 
        unsigned int min_cluster_size;

        // Scale parameter t. If present, a Gaussian smoothing based on the
        // scale value is executed before the clustering
        double scale;

        // Mean-shift bandwidth. Overrides the automatically calculated
        // bandwidth from scale when present.
        vector<T> ranges;

        // Contains a list of variables to exclude from scale-space filtering
        vector<std::string> exclude_from_scale_space_filtering;

        // The following kernel names are allowed:
        // 'uniform','epanechnikov','gauss','none'
        std::string kernel_name;

        // Name of the weight function to use. If adding your own
        // weight function, please give it a unique name and add
        // the code to create it to the class WeightFunctionFactory
        std::string weight_function_name;

        // Lower threshold for weight function filtering. Values at 
        // coordinates in Featurespace where the weight function is lower
        // than this value are omitted.
        T wwf_lower_threshold;

        // Upper threshold for weight function filtering. Values at 
        // coordinates in Featurespace where the weight function is higher
        // than this value are omitted.
        T wwf_upper_threshold;

        // If true, the weight function is written to the output file
        bool include_weight_in_result;

        // Important threshold determining how much cluster coverage
        // between previous and new clusters is required (in percent)
        // to qualify for merging/splitting. Leave alone if you don't
        // know what you are doing.
        T cluster_coverage_threshold;

        // Experimental flag. When true, the mean-shift calculation
        // is only performed for the spatial range of the Featurespace.
        bool spatial_range_only;

        // This string contains a string representation of the input parameters
        // after parsing for storing as attribute etc.
        string parameters;

        // When using previous clustering results, this contains the
        // filename of the previous cluster file. 
        std::string *previous_clusters_filename;

        // If present, the --previous-output file is used to consolidate 
        // current results. This is time consuming and has a propensity 
        // to form larger clusters. Use with discretion.
        bool postprocess_with_previous_output;

        // When convection filtering is used, this gives the index
        // of the variable to be filtered. Experimental. Ignore. 
        int convection_filter_index;

        // When this flag is true, neighboring clusters are coalesced 
        // according to their weight function response. Clusters with 
        // stronger weight function response 'absorb' neighboring clusters 
        // with lower response. Use cautious, very time consuming.
        bool coalesceWithStrongestNeighbour;

        // Verbosity of the processing chain. From 0 (silent) to 3 
        // (extremely verbose).
        Verbosity verbosity;

        // Flag indicating if the process was started with --inline-tracking
        // This has a number of consequences: cluster file will not be saved
        // after clustering (only after tracking) and the --previous-output
        // flag will be enforced
        bool inline_tracking;

        // ---------------------------------------------------------------
        // Replacement filtering.
        // ---------------------------------------------------------------

        // Variables indexed in this vector are subject to replacement
        // filtering. 
        vector<int> replacementFilterVariableIndex;

        // Replacement filter modes.
        map<int, typename ReplacementFilter<T>::ReplacementMode> replacementFilterModes;

        // Only used when mode is 'highest' or 'lowest'. Determines
        // the highest/lowest of the top/bottom percent are used.
        map<int, float> replacementFilterPercentages;

        // ---------------------------------------------------------------
        // CI-score weight function
        // ---------------------------------------------------------------
        // This is only relevant when working on the Oase2D composite
        // and CI score based tracking is wanted.

        // When CI-Score weight function is used, this file indicates 
        // the file to base the trend calculations on. Should be 15 minutes
        // in the past.
        std::string *ci_comparison_file;

        // Name of the protocluster file from 15 minutes ago, used for
        // trend calculation.
        std::string *ci_comparison_protocluster_file;

        // If this flag is true, the calculation is based solely on the
        // satellite info. Radar and Lightning are ignored.
        bool ci_satellite_only;

        // If this flag is true, the CI score calculation is performed
        // with the limits given in Mecikalski, John R., and Kristopher 
        // M. Bedka. “Forecasting Convective Initiation by Monitoring 
        // the Evolution of Moving Cumulus in Daytime GOES Imagery.” 
        // Monthly Weather Review 134, no. 1 (2006). doi:10.1175/MWR3062.1.
        bool ci_use_walker_mecikalski;

        // Scale of protoclusters
        T ci_protocluster_scale;

        // Minimum size of protoclusters
        int ci_protocluster_min_size;
        // -------------------------------------------------------------

#if WITH_VTK

        // Only present when VTK support is compiled in.
        // When true, a vtk file is written out containing the values
        // of the weight function. 
        bool write_weight_function;

        // Only present when VTK support is compiled in. 
        // When true, a vtk file containing the mean-shift vectors is
        // written out.
        bool write_meanshift_vectors;

        // Only present when VTK support is compiled in. 
        // When true, a vtk file containing the clusters and their weight
        // function values is written out.
        bool write_weight_response;

        // Only present when VTK support is compiled in. 
        // When true, a vtk file containing modes (mean-shift iteration end 
        // points) for all clusters is written out. Can be used to label
        // clusters.
        bool write_cluster_modes;

        // Only present when VTK support is compiled in. 
        // When true, a vtk file containing geometrical centers for all 
        // clusters is written out. Can be used to label clusters.
        bool write_cluster_centers;

        // Only present when VTK support is compiled in. 
        // When true, a vtk file containing the Featurespace after construction
        // and filtering.
        bool write_vtk;

        // Only present when VTK support is compiled in. All represented
        // variables in Featurespace are written out as .vtk files after
        // the processing.
        vector<std::string> vtk_variables;

#endif
    };

    /**
     * Data structure for holding context information during
     * a detection run.
     */
    template<typename T>
    struct detection_context_t
    {
        bool initialised;

        // TODO: to move away from NetCDF, this reference
        // should ultimately be removed
        NcFile *file;

        // This is what should be aimed for. An abstract data
        // store that could be anything.
        DataStore<T> *data_store;

        // Coordinate system
        CoordinateSystem<T> *coord_system;

        bool show_progress;
        timestamp_t timestamp;

        // Calculated automatically from scale or taken from the 
        // --range parameter if given
        vector<T> bandwidth;

        // Search parameters for mean-shift
        SearchParameters *search_params;

        // The main subject of interest
        FeatureSpace<T> *fs;

        // Automatically calculated 
        T kernel_width;
        T filter_width; // TODO: check if this is still needed?
        T decay; // Gauss filter decay threshold

        // Flag set to true if any of the weight function
        // filtering thresholds are set
        bool wwf_apply;

        // Scale-space filter. 
        ScaleSpaceFilter<T> *sf;

        // Weight function.
        WeightFunction<T> *weight_function;

        // Kernel for mean-shift 
        Kernel<T> *kernel;

        // K/D tree or other quick lookup index for mean-shift
        PointIndex<T> *index;

        // Resulting cluster list.
        typename ClusterList<T>::ptr clusters;

        // In case previous clusters are loaded, this contains those
        typename ClusterList<T>::ptr previous_clusters;
    };

    /** This class contains the tracking code.
    */
    template<typename T>
    class Detection
    {

    public:

        /** Constant representing 'no scale given. */
        static const double NO_SCALE;

        /**
         * Obtain a properly initialized default version of the
         * detection parameters. It is highly recommended to use
         * this method when creating detection parameters.
         * @return default detection parameters
         */
        static
        detection_params_t<T> defaultParams();

        /**
         * Initializes a context with defaults.
         * @param ctx
         */
        static
        void
        initialiseContext(detection_context_t<T> &ctx);

        /**
         * Initialize a detection context.
         * @param params
         * @return
         */
        static
        void
        initialiseContext(const detection_params_t<T> &params,
                          detection_context_t<T> &ctx);

        /**
         * Frees any memory allocated as a matter of parameter
         * parsing or processing.
         *
         * @param parameters
         * @param (initialised) context. If the context is not
         */
        static
        void
        cleanup(detection_params_t<T> &params,
                detection_context_t<T> &context);

        /**
         * Perform a detection run with the given parameters.
         * @param parameters
         * @param (initialised) context. If the context is not
         * initialized when the run starts, the method does it
         * for you.
         */
        static
        void
        run(const detection_params_t<T> &params,
            detection_context_t<T> &ctx);

    };
}

#endif	/* M3D_DETECTION_H */