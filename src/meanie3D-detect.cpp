

//
//  meanshift_clustering.cpp
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/exception/info.hpp>
#include <boost/smart_ptr.hpp>

#include <cf-algorithms/cf-algorithms.h>
#include <meanie3D/meanie3D.h>

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <exception>
#include <locale>
#include <limits>
#include <stdlib.h>
#include <netcdf>

using namespace std;
using namespace boost;
using namespace cfa::meanshift;
using namespace cfa::utils::visit;
using namespace cfa::utils::vectors;
using namespace netCDF;
using namespace m3D::utils::visit;
using namespace m3D::utils;
using namespace m3D::weights;

/** Feature-space data type */
typedef double FS_TYPE;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm,
        NcFile **filePtr,
        string &filename,
        string &output_filename,
        vector<NcDim> &dimensions,
        vector<NcVar> &dimension_variables,
        vector<NcVar> &variables,
        int &time_index,
        map<int, double> &lower_thresholds,
        map<int, double> &upper_thresholds,
        map<int, double> &replacement_values,
        double &scale,
        std::string &weight_function_name,
        FS_TYPE &wwf_lower_threshold,
        FS_TYPE &wwf_upper_threshold,
        int &convection_filter_index,
        string &parameters,
        vector<FS_TYPE> &ranges,
        std::string **previous_file,
        FS_TYPE &cluster_coverage_threshold,
        bool &coalesceWithStrongestNeighbour,
        bool &write_vtk,
        bool &write_weight_response,
        bool &write_cluster_centers,
        vector<size_t> &vtk_dimension_indexes,
        Verbosity &verbosity,
        unsigned int &min_cluster_size,
        vector<NcVar> &vtk_variables) {
    if (vm.count("file") == 0) {
        cerr << "Missing input file argument" << endl;

        exit(1);
    }

    filename = vm["file"].as<string>();

    try {
        output_filename = vm["output"].as<string>();
    }    catch (const boost::exception& e) {
        cerr << "Missing parameter -o " << endl;

        exit(-1);
    }

    // Open NetCDF file

    NcFile *file = NULL;

    try {
        file = new NcFile(filename, NcFile::read);
    }    catch (const netCDF::exceptions::NcException &e) {
        cerr << "ERROR opening file '" << filename << "' : " << e.what() << endl;
        exit(-1);
    }

    *filePtr = file;

    // Extract dimensions

    if (vm.count("dimensions") == 0) {
        cerr << "Missing parameter --dimensions" << endl;

        exit(1);
    }

    // parse dimension list

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    boost::char_separator<char> sep(",");

    tokenizer dim_tokens(vm["dimensions"].as<string>(), sep);

    for (tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter) {
        const char* name = (*tok_iter).c_str();

        dimensions.push_back(file->getDim(name));

        dimension_variables.push_back(file->getVar(name));
    }

    parameters = parameters + "dimensions=" + vm["dimensions"].as<string>() + " ";


    // parse variables

    if (vm.count("variables") == 0) {
        cerr << "Missing mandatory parameter --variables" << endl;

        exit(-1);
    }

    tokenizer var_tokens(vm["variables"].as<string>(), sep);

    for (tokenizer::iterator tok_iter = var_tokens.begin(); tok_iter != var_tokens.end(); ++tok_iter) {
        NcVar var = file->getVar(*tok_iter);

        if (var.isNull()) {
            cerr << "No variable '" << std::string(*tok_iter) << "' exists!" << endl;
            exit(-1);
        }

        variables.push_back(var);
    }

    parameters = parameters + "variables=" + vm["variables"].as<string>() + " ";

    // time 

    time_index = vm["time-index"].as<int>();

    // Convection filter index

    convection_filter_index = -1;

    if (vm.count("convection-filter-variable") > 0) {
        std::string cf_var_name = vm["convection-filter-variable"].as<string>();

        for (int index = 0; index < variables.size(); index++) {
            NcVar v = variables.at(index);

            if (v.getName() == cf_var_name) {
                // give the index in the complete feature-space vector
                // not just the list of variables handed in
                convection_filter_index = dimensions.size() + index;
            }
        }

        if (convection_filter_index < 0) {
            cerr << "Bad value for convection-filter-variable. Variable '" << cf_var_name << "' not a featurespace variable" << endl;
            exit(-1);
        }
    }


    // parse ranges if there

    if (vm.count("ranges") > 0) {
        tokenizer bw_tokens(vm["ranges"].as<string>(), sep);

        for (tokenizer::iterator tok_iter = bw_tokens.begin(); tok_iter != bw_tokens.end(); ++tok_iter) {
            const char* bw = (*tok_iter).c_str();

            ranges.push_back((FS_TYPE) strtod(bw, (char **) NULL));
        }

        if (ranges.size() != dimension_variables.size() + variables.size()) {
            cerr << "Please provide " << dimension_variables.size() + variables.size() << " bandwidth values" << endl;

            exit(1);
        }

        parameters = parameters + "ranges=" + vm["ranges"].as<string>();
    }

    // Lower Thresholds

    if (vm.count("lower-thresholds") > 0) {
        boost::char_separator<char> equals("=");

        tokenizer tokens(vm["lower-thresholds"].as<string>(), sep);

        for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            std::string pair = *tok_iter;

            tokenizer subtokens(pair, equals);

            tokenizer::iterator subtoken_iter = subtokens.begin();

            std::string variableName = *subtoken_iter;

            NcVar variable;

            for (size_t i = 0; i < variables.size(); i++) {
                if (variables[i].getName() == variableName) {
                    variable = variables[i];
                }
            }

            if (variable.isNull()) {
                cerr << "No variable named " << variableName << " found. Check --lower-thresholds parameter" << endl;

                exit(1);
            }

            subtoken_iter++;

            if (subtoken_iter == subtokens.end()) {
                cerr << "Missing threshold value for variable " << variableName << endl;

                exit(1);
            }

            const char* value = (*subtoken_iter).c_str();

            double doubleValue = strtod(value, (char **) NULL);

            lower_thresholds[variable.getId()] = doubleValue;
        }
    }

    // Upper Thresholds

    if (vm.count("upper-thresholds") > 0) {
        boost::char_separator<char> equals("=");

        tokenizer tokens(vm["upper-thresholds"].as<string>(), sep);

        for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            std::string pair = *tok_iter;

            tokenizer subtokens(pair, equals);

            tokenizer::iterator subtoken_iter = subtokens.begin();

            std::string variableName = *subtoken_iter;

            NcVar variable;

            for (size_t i = 0; i < variables.size(); i++) {
                if (variables[i].getName() == variableName) {
                    variable = variables[i];
                }
            }

            if (variable.isNull()) {
                cerr << "No variable named " << variableName << " found. Check --upper-thresholds parameter" << endl;

                exit(1);
            }

            subtoken_iter++;

            if (subtoken_iter == subtokens.end()) {
                cerr << "Missing threshold value for variable " << variableName << endl;
                exit(1);
            }

            const char* value = (*subtoken_iter).c_str();

            double doubleValue = strtod(value, (char **) NULL);

            upper_thresholds[variable.getId()] = doubleValue;
        }
    }

    // Fill values

    if (vm.count("replacement-values") > 0) {
        boost::char_separator<char> equals("=");

        tokenizer tokens(vm["replacement-values"].as<string>(), sep);

        for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
            std::string pair = *tok_iter;

            tokenizer subtokens(pair, equals);

            tokenizer::iterator subtoken_iter = subtokens.begin();

            std::string variableName = *subtoken_iter;

            NcVar variable;

            for (size_t i = 0; i < variables.size(); i++) {
                if (variables[i].getName() == variableName) {
                    variable = variables[i];
                }
            }

            if (variable.isNull()) {
                cerr << "No variable named " << variableName << " found. Check --replacement-values parameter" << endl;

                exit(1);
            }

            subtoken_iter++;

            if (subtoken_iter == subtokens.end()) {
                cerr << "Missing replacement value for variable " << variableName << endl;
                exit(1);
            }

            const char* value = (*subtoken_iter).c_str();

            double doubleValue = strtod(value, (char **) NULL);

            replacement_values[variable.getId()] = doubleValue;
        }
    }


    // Scale parameter

    scale = vm["scale"].as<double>();

    // Weight Function

    weight_function_name = vm["weight-function-name"].as<string>();

    if (!(weight_function_name == "default" || weight_function_name == "inverse" || weight_function_name == "oase" || weight_function_name == "pow10")) {
        cerr << "Illegal weight function name " << weight_function_name << ". Only 'default','inverse','pow10' or 'oase' are known." << endl;
        exit(1);
    }

    wwf_lower_threshold = vm["wwf-lower-threshold"].as<FS_TYPE>();

    wwf_upper_threshold = vm["wwf-upper-threshold"].as<FS_TYPE>();

    // Coalescence?

    coalesceWithStrongestNeighbour = vm.count("coalesce-with-strongest-neighbour") > 0;

    // VTK output?

    write_vtk = vm.count("write-clusters-as-vtk") > 0;

    write_cluster_centers = vm.count("write-cluster-centers") > 0;

    write_weight_response = vm.count("write-cluster-weight-response") > 0;

    // VTK dimension mapping

    if (vm.count("vtk-dimensions") > 0) {
        // parse dimension list

        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        boost::char_separator<char> sep(",");

        string str_value = vm["vtk-dimensions"].as<string>();

        tokenizer dim_tokens(str_value, sep);

        for (tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter) {
            const char* name = (*tok_iter).c_str();

            NcDim dim = file->getDim(name);

            vector<NcDim>::const_iterator fi = find(dimensions.begin(), dimensions.end(), dim);

            if (fi == dimensions.end()) {
                cerr << "--vtk-dimension parameter " << dim.getName() << " is not part of --dimensions" << endl;
                exit(-1);
            }

            size_t index = fi - dimensions.begin();

            vtk_dimension_indexes.push_back(index);
        }

        if (vtk_dimension_indexes.size() != dimensions.size()) {
            cerr << "The number of vtk-dimensions must be identical to dimensions" << endl;

            exit(-1);
        }
    }

    // Previous file

    if (vm.count("previous-file") > 0) {
        std::string previous = vm["previous-file"].as<string>();

        boost::filesystem::path previous_path(previous);

        if (boost::filesystem::exists(previous_path) && boost::filesystem::is_regular_file(previous_path)) {
            *previous_file = new std::string(previous);
        } else {
            cerr << "Illegal value for parameter --previous-file: does not exist or is no regular file" << endl;
        }
    }

    cluster_coverage_threshold = vm["previous-cluster-coverage-threshold"].as<FS_TYPE>();

    // Verbosity

    unsigned short vb = vm["verbosity"].as<unsigned short>();

    if (vb > VerbosityAll) {
        cerr << "Illegal value for parameter --verbosity. Only values from 0 .. 3 are allowed" << endl;

        exit(-1);
    } else {
        verbosity = (Verbosity) vb;
    }

    min_cluster_size = vm["min-cluster-size"].as<unsigned int>();

    // vtk-variables

    if (vm.count("write-variables-as-vtk") > 0) {
        tokenizer bw_tokens(vm["write-variables-as-vtk"].as<string>(), sep);

        for (tokenizer::iterator tok_iter = bw_tokens.begin(); tok_iter != bw_tokens.end(); ++tok_iter) {
            const char* bw = (*tok_iter).c_str();

            try {
                NcVar var = file->getVar(bw);

                if (var.isNull()) {
                    cerr << "Can't open variable " << bw << " from NetCDF file. Check --write-variables-as-vtk" << endl;
                    exit(-1);
                }

                vtk_variables.push_back(var);
            }            catch (const netCDF::exceptions::NcException &e) {
                cerr << "Can't find variable " << bw << " from NetCDF file. Check --write-variables-as-vtk" << endl;
                exit(-1);
            }
        }
    }
}

/**
 *
 *
 */
int main(int argc, char** argv) {
    using namespace cfa::utils::timer;
    using namespace cfa::utils::visit;
    using namespace m3D;

    // Declare the supported options.

    program_options::options_description desc("Options");
    desc.add_options()
            ("help,h", "produce help message")
            ("file,f", program_options::value<string>(), "CF-Metadata compliant NetCDF-file")
            ("output,o", program_options::value<string>(), "Name of output file for clustering results")
            ("dimensions,d", program_options::value<string>(), "Comma-separatred list of the dimensions to be used. The program expects dimension variables with identical names.")
            ("vtk-dimensions", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
            ("variables,v", program_options::value<string>(), "Comma-separated variables used to construct feature space. Do not include dimension variables")
            ("time-index,t", program_options::value<int>()->default_value(-1), "Index of the point in time you wish to use in files with a time dimension. -1 means no time dimension in the variables.")
            ("lower-thresholds", program_options::value<string>(), "Comma-separated list var1=val,var2=val,... of lower tresholds. Values below this are ignored when constructing feature space")
            ("upper-thresholds", program_options::value<string>(), "Comma-separated list var1=val,var2=val,... of lower tresholds. Values above this are ignored when constructing feature space")
            ("replacement-values", program_options::value<string>(), "Comma-separated list var1=val,var2=val,... of values to replace missing values with in feature space construction. If no replacement value is specified while even one variable is out of valid range at one point, the whole point is discarded")
            ("weight-function-name,w", program_options::value<string>()->default_value("default"), "default,inverse,pow10 or oase")
            ("wwf-lower-threshold", program_options::value<FS_TYPE>()->default_value(0.05), "Lower threshold for weight function filter. Defaults to 0.05 (5%)")
            ("wwf-upper-threshold", program_options::value<FS_TYPE>()->default_value(std::numeric_limits<FS_TYPE>::max()), "Upper threshold for weight function filter. Defaults to std::numeric_limits::max()")
            ("coalesce-with-strongest-neighbour", "Clusters are post-processed, coalescing each cluster with their strongest neighbour")
            ("convection-filter-variable,c", program_options::value<string>(), "Name the variable to eliminate all but the points marked as convective according to the classification scheme by Steiner,Houza & Yates (1995) in.")
            ("scale,s", program_options::value<double>()->default_value(NO_SCALE), "Scale parameter to pre-smooth the data with. Filter size is calculated from this automatically.")
            ("filter-size,l", program_options::value<double>()->default_value(0.0), "Scale parameter to pre-smooth the data with. Scale parameter is calculated from this automatically.")
            ("ranges,r", program_options::value<string>(), "Override the automatic bandwidth calculation with a set of given bandwidths. Use in the order of (dim1,...dimN,var1,...,varN).")
            ("min-cluster-size,m", program_options::value<unsigned int>()->default_value(1u), "Discard clusters smaller than this number of points.")
            ("previous-file,p", program_options::value<string>(), "Optional file containing the clustering results from the previous timeslice. Helps to keep the clustering more stable over time.")
            ("previous-cluster-coverage-threshold", program_options::value<double>()->default_value(0.66), "Minimum overlap in percent between current and previous clusters to be taken into consideration. Defaults to 2/3 (0.66)")
            ("write-variables-as-vtk", program_options::value<string>(), "Comma separated list of variables that should be written out as VTK files (after applying scale/threshold)")
            ("write-clusters-as-vtk", "write clusters out in .vtk file format additionally (useful for visualization with visit for example)")
            ("write-cluster-centers", "write cluster centers out in .vtk file format (useful for visualization with visit for example)")
            ("write-cluster-weight-response", "write out the clusters with weight responses as value")
            ("verbosity", program_options::value<unsigned short>()->default_value(1), "Verbosity level [0..3], 0=silent, 1=normal, 2=show details, 3=show all details). Default is 1.")
            ;

    program_options::variables_map vm;
    try {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    }    catch (std::exception &e) {
        cerr << "Error parsing command line: " << e.what() << endl;
        cerr << "Check meanie3D-detect --help for command line options" << endl;
        exit(-1);
    }

    if (vm.count("help") == 1 || argc < 2) {
        cout << desc << "\n";
        return 1;
    }

    // Evaluate user input

    NcFile *file = NULL;
    string filename;
    string output_filename;
    string parameters;
    vector<NcDim> dimensions;
    vector<size_t> vtk_dimension_indexes;
    vector<NcVar> dimension_variables;
    vector<NcVar> variables;
    int time_index = 0;
    vector<double> ranges;
    unsigned int min_cluster_size = 1;
    map<int, double> lower_thresholds; // ncvar.id / value
    map<int, double> upper_thresholds; // ncvar.id / value
    map<int, double> replacement_values; // ncvar.id / value
    FS_TYPE wwf_lower_threshold;
    FS_TYPE wwf_upper_threshold;
    std::string weight_function_name = "default";
    std::string *previous_file = NULL;
    FS_TYPE cluster_coverage_threshold = 0.66;
    int convection_filter_index = -1;
    bool coalesceWithStrongestNeighbour = false;
    double scale = NO_SCALE;

    vector<NcVar> vtk_variables;
    SearchParameters *search_params = NULL;
    bool write_vtk = false;
    bool write_weight_response = false;
    bool write_cluster_centers = false;
    Verbosity verbosity = VerbosityNormal;

    timestamp_t timestamp;

    FS_TYPE kernel_width = 0.0;

    try {
        parse_commmandline(vm,
                &file,
                filename,
                output_filename,
                dimensions,
                dimension_variables,
                variables,
                time_index,
                lower_thresholds,
                upper_thresholds,
                replacement_values,
                scale,
                weight_function_name,
                wwf_lower_threshold,
                wwf_upper_threshold,
                convection_filter_index,
                parameters,
                ranges,
                &previous_file,
                cluster_coverage_threshold,
                coalesceWithStrongestNeighbour,
                write_vtk,
                write_weight_response,
                write_cluster_centers,
                vtk_dimension_indexes,
                verbosity,
                min_cluster_size,
                vtk_variables);

        // Make the mapping known to the visualization routines

        ::cfa::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        ::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;

        // Select the correct point factory
        PointFactory<FS_TYPE>::set_instance(new M3DPointFactory<FS_TYPE>());

        // Get timestamp

        // Tracking comparison (solve that generically sometime)
        if (boost::contains(filename, "rico.out.xy.")) {
            // time for this is days since simulation start
            double time_in_days = get_time<double>(filename, time_index);

            // a day has 24 * 60 * 60 seconds
            timestamp = (long) round(time_in_days * 24.0 * 60.0 * 60.0);
        } else {
            // Seconds since epoch
            timestamp = get_time<long>(filename, time_index);
        }
    }    catch (const std::exception &e) {
        cerr << e.what() << endl;
        exit(-1);
    }

    bool show_progress = (verbosity > VerbositySilent);

    if (verbosity > VerbositySilent) {
        cout << "----------------------------------------------------" << endl;
        cout << "Meanie3D-detect" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << endl;

        cout << "Command line options:" << endl;

        cout << "\tinput file = " << filename << endl;

        if (time_index >= 0) {
            cout << "\tusing point in time at index " << time_index << " (timestamp=" << timestamp << ")" << endl;
        } else {
            cout << "\ttime is not a variable dimension. Using timestamp " << timestamp << endl;
        }

        cout << "\tdimensions = " << vm["dimensions"].as<string>() << endl;

        cout << "\tvariables = " << vm["variables"].as<string>() << endl;

        if (!ranges.empty()) {
            cout << "\tranges = " << ranges << endl;

            // calculate kernel width as average of ranges

            kernel_width = boost::numeric_cast<FS_TYPE>(0.0);

            for (size_t i = 0; i < dimensions.size(); i++) {
                kernel_width += ranges[i];
            }

            kernel_width = kernel_width / boost::numeric_cast<FS_TYPE>(dimensions.size());
        } else {
            cout << "\tautomatic bandwidth selection" << endl;
        }

        if (!lower_thresholds.empty()) {
            cout << "\tusing lower thresholds " << vm["lower-thresholds"].as<string>() << endl;
        }

        if (!upper_thresholds.empty()) {
            cout << "\tusing upper thresholds " << vm["upper-thresholds"].as<string>() << endl;
        }

        if (scale != NO_SCALE) {
            double width = sqrt(ceil(-2.0 * scale * log(0.01))) / 2;

            if (ranges.empty()) {
                kernel_width = width;
            }

            cout << "\tpre-smoothing data with scale parameter " << scale << " (kernel width = " << width << ")" << endl;
        } else {
            cout << "\tno scale-space smoothing" << endl;
        }

        cout << "\tweight-function:" << weight_function_name << endl;
        cout << "\t\tlower weight-function threshold: " << wwf_lower_threshold << endl;
        cout << "\t\tupper weight-function threshold: " << wwf_upper_threshold << endl;

        if (previous_file != NULL)
        {
            cout << "\tprevious file:" << *previous_file << endl;

            if (vm.count("previous-cluster-coverage-threshold") > 0) {
                cout << "\tprevious file coverage threshold:" << cluster_coverage_threshold << endl;
            }
        }

        if (min_cluster_size > 1) {
            cout << "\tminimum cluster size = " << min_cluster_size << endl;
        }
        
        cout << "\tcoalesce results with strongest neighbor: " << (coalesceWithStrongestNeighbour ? "yes":"no") << endl;

        cout << "\toutput written to file: " << output_filename << endl;

        if (!vtk_variables.empty()) {
            cout << "\twriting out these variables as vtk after processing:" << vm["write-variables-as-vtk"].as<string>() << endl;
        }

        cout << "\tclusters written as vtk: " << (write_vtk ? "yes" : "no") << endl;
        cout << "\twrite cluster centers as vtk: " << (write_cluster_centers ? "yes" : "no") << endl;
        cout << "\twriting cluster weights as vtk: " << (write_weight_response ? "yes" : "no") << endl;

        if (verbosity > VerbosityNormal) {
            cout << endl;
            cout << "Compiled options (use -D to switch them on and off during compile time)" << endl;

#if WRITE_BANDWIDTH
            cout << "\tWRITE_BANDWIDTH=1" << endl;
#else
            cout << "\tWRITE_BANDWIDTH=0" << endl;
#endif

#if WRITE_CLUSTER_MEANSHIFT
            cout << "\tWRITE_CLUSTER_MEANSHIFT=1" << endl;
#else
            cout << "\tWRITE_CLUSTER_MEANSHIFT=0" << endl;
#endif

#if WRITE_CLUSTER_MODES
            cout << "\tWRITE_CLUSTER_MODES=1" << endl;
#else
            cout << "\tWRITE_CLUSTER_MODES=0" << endl;
#endif

#if WRITE_FEATURESPACE
            cout << "\tWRITE_FEATURESPACE=1" << endl;
#else
            cout << "\tWRITE_FEATURESPACE=0" << endl;
#endif

#if WRITE_INDEX
            cout << "\tWRITE_INDEX=1" << endl;
#else
            cout << "\tWRITE_INDEX=0" << endl;
#endif

#if WRITE_ITERATION_ORIGINS
            cout << "\tWRITE_ITERATION_ORIGINS=1" << endl;
#else
            cout << "\tWRITE_ITERATION_ORIGINS=0" << endl;
#endif

#if WRITE_MEANSHIFT_SAMPLES
            cout << "\tWRITE_MEANSHIFT_SAMPLES=1" << endl;
#else
            cout << "\tWRITE_MEANSHIFT_SAMPLES=0" << endl;
#endif

#if WRITE_MEANSHIFT_VECTORS
            cout << "\tWRITE_MEANSHIFT_VECTORS=1" << endl;
#else
            cout << "\tWRITE_MEANSHIFT_VECTORS=0" << endl;
#endif

#if WRITE_MEANSHIFT_WEIGHTS
            cout << "\tWRITE_MEANSHIFT_WEIGHTS=1" << endl;
#else
            cout << "\tWRITE_MEANSHIFT_WEIGHTS=0" << endl;
#endif

#if WRITE_MODES
            cout << "\tWRITE_MODES=1" << endl;
#else
            cout << "\tWRITE_MODES=0" << endl;
#endif

#if WRITE_OFF_LIMITS_MASK
            cout << "\tWRITE_OFF_LIMITS_MASK=1" << endl;
#else
            cout << "\tWRITE_OFF_LIMITS_MASK=0" << endl;
#endif

#if WRITE_WEIGHT_FUNCTION
            cout << "\tWRITE_WEIGHT_FUNCTION=1" << endl;
#else
            cout << "\tWRITE_WEIGHT_FUNCTION=0" << endl;
#endif

#if WRITE_ZEROSHIFT_CLUSTERS
            cout << "\tWRITE_ZEROSHIFT_CLUSTERS=1" << endl;
#else
            cout << "\tWRITE_ZEROSHIFT_CLUSTERS=0" << endl;
#endif
            cout << endl;
        }
    }

    // Construct Featurespace

    // Coordinate system

    CoordinateSystem<FS_TYPE> *coord_system = new CoordinateSystem<FS_TYPE>(dimensions, dimension_variables);

    // Feature Space

    start_timer();

    // TODO: threshold should be applied as a filter somehow
    // this approach seems a little half-cocked

    FeatureSpace<FS_TYPE> *fs = new FeatureSpace<FS_TYPE>(filename,
            coord_system,
            variables,
            time_index,
            lower_thresholds,
            upper_thresholds,
            replacement_values,
            show_progress);

    // if ranges are not explicitly given, figure the out
    // based on scale

    if (ranges.empty()) {
        // spatial range first

        FS_TYPE t = (scale == NO_SCALE) ? 1.0 : scale;

        FS_TYPE filter_width = sqrt(ceil(-2.0 * t * log(0.01))) / 2.0;

        for (size_t i = 0; i < dimensions.size(); i++)
        {
            ranges.push_back(filter_width);
        }

        // value range

        for (size_t i = 0; i < variables.size(); i++) {
            FS_TYPE range = fs->valid_max().at(i) - fs->valid_min().at(i);

            ranges.push_back(range);
        }

        cout << "Automatically calculated bandwidth:" << ranges << endl;
    }

    search_params = new RangeSearchParams<FS_TYPE>(ranges);


#if WRITE_OFF_LIMITS_MASK
    fs->off_limits()->write("off_limits.vtk", "off_limits");
#endif

    // used in writing out debug data
    boost::filesystem::path path(filename);

    // Convection Filter?

    if (convection_filter_index >= 0) {
        ConvectionFilter<FS_TYPE> convection_filter(ranges, convection_filter_index, show_progress);
        convection_filter.apply(fs);
    }

    // Scale-Space Smoothing

    WeightFunction<FS_TYPE> *weight_function = NULL;

    // Scale-Space smoothing

    if (scale != NO_SCALE) {
        // TODO: make decay a parameter or at least a constant

        FS_TYPE decay = 0.01;

        ScaleSpaceFilter<FS_TYPE> sf(scale, fs->coordinate_system->resolution(), decay, show_progress);

        sf.apply(fs);

        if (verbosity > VerbositySilent)
            cout << endl << "Constructing " << weight_function_name << " weight function ...";

        if (weight_function_name == "oase")
        {
            weight_function = new OASEWeightFunction<FS_TYPE>(fs, sf.get_filtered_min(), sf.get_filtered_max());
        }
        else if (weight_function_name == "inverse")
        {
            weight_function = new InverseDefaultWeightFunction<FS_TYPE>(fs, sf.get_filtered_min(), sf.get_filtered_max());
        }
        else if (weight_function_name == "pow10")
        {
            weight_function = new EXP10WeightFunction<FS_TYPE>(fs);
        }
        else
        {
            weight_function = new DefaultWeightFunction<FS_TYPE>(fs, sf.get_filtered_min(), sf.get_filtered_max());
        }

#if WRITE_WEIGHT_FUNCTION
        std::string wfname = path.filename().stem().string() + "-weights";
        ::cfa::utils::VisitUtils<FS_TYPE>::write_weight_function_response(wfname, fs, weight_function);
#endif
        if (verbosity > VerbositySilent)
            cout << " done." << endl;

        // Apply weight function filter
        // TODO: make that a command line parameter

        WeightThresholdFilter<FS_TYPE> wtf(weight_function, wwf_lower_threshold, wwf_upper_threshold, true);

        wtf.apply(fs);

        if (verbosity > VerbositySilent)
            cout << "Filtered featurespace contains " << fs->count_original_points() << " original points " << endl;

    }
    else
    {
        if (verbosity > VerbositySilent)
            cout << endl << "Constructing " << weight_function_name << " weight function ...";

        if (weight_function_name == "oase")
        {
            weight_function = new OASEWeightFunction<FS_TYPE>(fs);
        }
        else if (weight_function_name == "inverse")
        {
            weight_function = new InverseDefaultWeightFunction<FS_TYPE>(fs, lower_thresholds, upper_thresholds);
        }
        else if (weight_function_name == "pow10")
        {
            weight_function = new EXP10WeightFunction<FS_TYPE>(fs);
        }
        else
        {
            weight_function = new DefaultWeightFunction<FS_TYPE>(fs);
        }

#if WRITE_WEIGHT_FUNCTION
        std::string wfname = path.filename().stem().string() + "-weights";
        ::cfa::utils::VisitUtils<FS_TYPE>::write_weight_function_response(wfname, fs, weight_function);
#endif

        // Apply weight function filter

        WeightThresholdFilter<FS_TYPE> wtf(weight_function, wwf_lower_threshold, wwf_upper_threshold, true);

        wtf.apply(fs);

        if (verbosity > VerbositySilent)
            cout << "Filtered featurespace contains " << fs->count_original_points() << " original points " << endl;

        if (verbosity > VerbositySilent)
            cout << " done." << endl;
    }

    if (!vtk_variables.empty()) {
        string filename_only = boost::filesystem::path(filename).filename().string();

        boost::filesystem::path destination_path = boost::filesystem::path(".");

        destination_path /= filename_only;

        //destination_path.replace_extension("vtk");

        destination_path.replace_extension();

        string dest_path = destination_path.generic_string();

        if (verbosity > VerbositySilent)
            cout << "Writing featurespace-variables ...";

        cfa::utils::VisitUtils<FS_TYPE>::write_featurespace_variables_vtk(dest_path, fs, vtk_variables, false);

        if (verbosity > VerbositySilent)
            cout << " done." << endl;
    }

    if (verbosity == VerbosityAll)
        fs->print();

    // Calculate termcrit_epsilon if necessary

#if WRITE_MEANSHIFT_WEIGHTS
    vector<FS_TYPE> sample_point(2);
    sample_point[0] = -43.4621658325195;
    sample_point[1] = -4383.64453125;
    fs->weight_sample_points.push_back(sample_point);
#endif

    PointIndex<FS_TYPE> *index = PointIndex<FS_TYPE>::create(fs);

    //PointIndex<FS_TYPE> *index = PointIndex<FS_TYPE>::create(fs,PointIndex<FS_TYPE>::IndexTypeRectilinearGrid);

    //
    // Simple clustering
    //

    ClusterOperation<FS_TYPE> cop(fs, index);

    // estimate kernel size

    Kernel<FS_TYPE> *kernel = new UniformKernel<FS_TYPE>(kernel_width);

    // Kernel<FS_TYPE> *kernel = new EpanechnikovKernel<FS_TYPE>(1.0);

    ClusterList<FS_TYPE> clusters = cop.cluster(search_params, kernel, weight_function, coalesceWithStrongestNeighbour, show_progress);

    // Sanity check

    // clusters.sanity_check( fs );

    // Number the result sequentially to make it easier to follow
    // previous results in comparison

    clusters.retag_identifiers();

    // Axe weenies

    clusters.apply_size_threshold(min_cluster_size);

    if (verbosity == VerbosityAll) {
        clusters.print();
    }

    // Collate with previous clusters, if those are provided
    
    if (previous_file != NULL)
    {
        cout << endl << "Collating with previous results. Have " << clusters.clusters.size() << " clusters:" << endl;
        clusters.print();

        try {
            ClusterList<FS_TYPE>::ptr previous = ClusterList<FS_TYPE>::read(*previous_file);

            ClusterUtils<FS_TYPE> cluster_filter(cluster_coverage_threshold);

            cluster_filter.filter_with_previous_clusters(previous, &clusters, weight_function, verbosity);
        }        catch (const std::exception &e) {
            cerr << "ERROR reading previous cluster file: " << e.what() << endl;
            exit(-1);
        }
        
        cout << endl << "Done. Have " << clusters.clusters.size() << " clusters:" << endl;
        clusters.print();

    }

    // Announce final results

    if (verbosity > VerbositySilent) {
        cout << endl << "Final result: found " << clusters.clusters.size() << " objects: " << endl;

        clusters.print();
    }

    // Write out the cluster list

    if (write_vtk && clusters.clusters.size() > 0) {
        //::m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtk(clusters, coord_system, path.filename().string());
        ::m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtr(&clusters, coord_system, path.filename().string());
    }

#if WRITE_CLUSTER_MODES
    // MODES are needed for tagging with IDs
    string modes_path = path.filename().stem().string() + "-clusters_modes.vtk";
    ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_modes_vtk(modes_path, clusters.clusters, true);
#endif

    if (write_cluster_centers) {
        string centers_path = path.filename().stem().string() + "-clusters_centers.vtk";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_geometrical_cluster_centers_vtk(centers_path, clusters.clusters);
    }

#if WRITE_CLUSTER_MEANSHIFT
    string shifts_path = path.filename().stem().string() + "-clusters_shifts";
    ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_meanshift_vtk(shifts_path, clusters.clusters);
#endif

    if (write_weight_response && clusters.clusters.size() > 0) {
        string wr_path = path.filename().stem().string() + "-clusters_weight";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_weight_response_vtk(wr_path, clusters.clusters, weight_function, false);
    }

    if (verbosity > VerbositySilent)
        cout << "Writing clusters to NetCDF file " << output_filename << " ..." << endl;

    // Before writing, set the timestamp!!

    clusters.timestamp = timestamp;

    clusters.write(output_filename);

    if (verbosity > VerbositySilent)
        cout << "done." << endl;

    // mop up

    delete kernel;

    delete index;

    delete coord_system;

    delete fs;

    return 0;
}




