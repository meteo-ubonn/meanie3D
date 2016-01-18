//
//  meanie3D-detect
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
#include <string>
#include <netcdf>
#include <meanie3D/tracking/tracking.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;
using namespace m3D::utils;

/** Feature-space data type */
typedef double FS_TYPE;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm, tracking_param_t &params)
{
    if (vm.count("version") != 0) {
        cout << m3D::VERSION << endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("previous") == 0) {
        cerr << "Missing previous cluster file argument --previous" << endl;
        exit(1);
    }

    try {
        params.previous_filename = vm["previous"].as<string>();
        NcFile previous_cluster_file(params.previous_filename, NcFile::read);
    } catch (const netCDF::exceptions::NcException &e) {
        cerr << "Exception opening file " << params.previous_filename << ":" << endl;
        cerr << "FATAL:" << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    try {
        params.current_filename = vm["current"].as<string>();
        NcFile current_cluster_file(params.current_filename, NcFile::write);
    } catch (const netCDF::exceptions::NcException &e) {
        cerr << "FATAL:exception opening file " << params.current_filename << ":" << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Verbosity
    unsigned short vb = vm["verbosity"].as<unsigned short>();
    if (vb > VerbosityAll) {
        cerr << "Illegal value for parameter --verbosity. Only values from 0 .. 3 are allowed" << endl;
        exit(EXIT_FAILURE);
    } else {
        params.verbosity = (Verbosity) vb;
    }

    // max-speed
    FS_TYPE speed = vm["max-speed"].as<FS_TYPE>();
    params.maxVelocity= ::units::values::meters_per_second(speed);

    // max time
    FS_TYPE time = vm["max-time"].as<FS_TYPE>();
    params.max_deltaT= ::units::values::s(time);

    // continue ID?
    params.continueIDs = !(vm.count("discontinue-id-in-merge-and-split") > 0);
    params.useDisplacementVectors = vm.count("use-displacement-vectors") > 0;

    // merge/split continuation threshold
    params.mergeSplitContinuationThreshold = vm["merge-split-continuation-threshold"].as<FS_TYPE>();

    // merge/split threshold
    params.mergeSplitThreshold = vm["merge-split-threshold"].as<FS_TYPE>();

    // Weights
    params.range_weight = vm["wr"].as<FS_TYPE>();
    params.size_weight = vm["ws"].as<FS_TYPE>();
    params.correlation_weight = vm["wt"].as<FS_TYPE>();

    // tracking variable
    if (vm.count("tracking-variable") > 0) {
        params.tracking_variable = vm["tracking-variable"].as<string>();
    }

#if WITH_VTK
    // --write-vtk
    params.write_vtk = vm.count("write-vtk") > 0;
    // VTK dimension mapping
    // TODO: This code should be made generic and moved into utils.
    if (vm.count("vtk-dimensions") > 0) {
        // Open the file for reading once more
        NcFile *file = new NcFile(params.current_filename, NcFile::read);
        // Read "featurespace_dimensions"
        string fs_dimensions;
        file->getAtt("featurespace_dimensions").getValues(fs_dimensions);
        vector<string> fs_dim_names = vectors::from_string<string>(fs_dimensions);
        // parse dimension list
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep(",");
        string str_value = vm["vtk-dimensions"].as<string>();
        tokenizer dim_tokens(str_value, sep);
        for (tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter) {
            string name = *tok_iter;
            // Boost::tokenizer has a horrible bug that only hits when
            // Optimization is above -O2. It appends brackets to the
            // beginning or end of the token. Totally fucked up shit
            bool found_it = false;
            for (size_t pi = 0; pi < fs_dim_names.size(); pi++) {
                std::string dim_name = fs_dim_names[pi];
                if (name.compare(dim_name) == 0) {
                    params.vtk_dimension_indexes.push_back((size_t) pi);
                    found_it = true;
                    break;
                }
            }
            if (!found_it) {
                cerr << "FATAL:invalid dimension '" << name << "'. Check parameter --vtk-dimensions" << endl;
                exit(EXIT_FAILURE);
            }
        }
        delete file;
        VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = params.vtk_dimension_indexes;
    }

#endif
}

/**
 *
 *
 */
int main(int argc, char** argv) {
    // Declare the supported options.

    tracking_param_t params = Tracking<FS_TYPE >::defaultParams();

    program_options::options_description desc("Options");
    desc.add_options()
            ("help,h", "produce help message")
            ("version", "print version information and exit")
            ("previous,p", program_options::value<string>(), "Previous cluster file (netCDF)")
            ("current,c", program_options::value<string>(), "Current cluster file (netCDF)")
            ("tracking-variable,t", program_options::value<string>()->default_value(params.tracking_variable), "Variable used for histogram correlation. Must be specified when histogram weight --wt is not zero")
            ("wr", program_options::value<FS_TYPE>()->default_value(params.range_weight), "Weight for range correlation [0..1]")
            ("ws", program_options::value<FS_TYPE>()->default_value(params.size_weight), "Weight for size correlation [0..1]")
            ("wt", program_options::value<FS_TYPE>()->default_value(params.correlation_weight), "Weight for histogram rank correlation [0..1]")
            ("use-displacement-vectors,v", "If present, the algorithm uses the displacement vectors from the previous tracking result (if present) to shift clusters from the previous file to improve tracking (experimental).")
            ("merge-split-threshold", program_options::value<FS_TYPE>()->default_value(params.mergeSplitThreshold), "Percentage of area covered between previous/new clusters for split/merge calculation")
            ("merge-split-continuation-threshold", program_options::value<FS_TYPE>()->default_value(params.mergeSplitContinuationThreshold), "Minimum percentage of area covered between previous/new clusters to continue ID")
            ("discontinue-id-in-merge-and-split,d", "If present, the tracking discontinues cluster IDs when merging and splitting. Otherwise the largest candidate carries the ID on if the overlap exceeds --merge-split-continuation-threshold.")
            ("max-speed", program_options::value<FS_TYPE>()->default_value(params.maxVelocity.get()), "Maximum allowed object speed (m/s)")
            ("max-time", program_options::value<FS_TYPE>()->default_value(params.max_deltaT.get()), "Maximum allowed time difference between files (seconds)")
            ("max-size-deviation", program_options::value<double>()->default_value(params.max_size_deviation), "Maximum allowed difference in sizes (number of points) between clusters from previous and current file in percent [0..1]")
#if WITH_VTK
            ("write-vtk,k", "Write out the clusters as .vtk files for visit")
            ("vtk-dimensions", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
#endif
            ("verbosity", program_options::value<unsigned short>()->default_value(1), "Verbosity level [0..3], 0=silent, 1=normal, 2=show details, 3=show all details). Default is 1.")
            ;

    program_options::variables_map vm;

    try {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    } catch (std::exception &e) {
        cerr << "FATAL:error parsing command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help") == 1 || argc < 2) {
        cout << desc << "\n";
        return 1;
    }

    // Evaluate user input
    try {
        parse_commmandline(vm, params);
    } catch (const std::exception &e) {
        cerr << "FATAL:" << e.what() << endl;
        exit(EXIT_FAILURE);
    }


    if (params.verbosity > VerbositySilent) {
        cout << "----------------------------------------------------" << endl;
        cout << "meanie3D-track" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << endl;
        cout << "Command line options:" << endl;
        cout << "\tprevious cluster file = " << params.previous_filename << endl;
        cout << "\tcurrent cluster file = " << params.current_filename << endl;
        cout << "\tvariable name for histogram comparison: " << params.tracking_variable << endl;
        cout << "\tcorrelation weights: wr=" << params.range_weight << " ws=" << params.size_weight << " wt=" << params.correlation_weight << endl;
        cout << "\tmaximum speed: " << params.maxVelocity << " [m/s]" << endl;
        cout << "\tmaximum time difference: " << params.max_deltaT << " [seconds]" << endl;
        cout << "\tmerge/split threshold: " << params.mergeSplitThreshold << endl;
        cout << "\tcontinue id in merge/split: " << (params.continueIDs ? "yes" : "no") << endl;
        cout << "\tmerge/split id continuation threshold: " << params.mergeSplitContinuationThreshold << endl;
        cout << "\tusing displacement vectors to shift previous clusters: " << (params.useDisplacementVectors ? "yes" : "no") << endl;
#if WITH_VTK
        cout << "\twriting results out as vtk:" << (params.write_vtk ? "yes" : "no") << endl;
#endif
        cout << endl;
    }

    // Read previous clusters
    if (params.verbosity >= VerbosityNormal) start_timer("Reading " + params.previous_filename+ " ... ");
    ClusterList<FS_TYPE>::ptr previous = ClusterList<FS_TYPE>::read(params.previous_filename);
    if (params.verbosity >= VerbosityNormal) stop_timer("done");

    // Read current clusters
    CoordinateSystem<FS_TYPE> *cs;
    if (params.verbosity >= VerbosityNormal) start_timer("Reading " + params.current_filename + " ... ");
    ClusterList<FS_TYPE>::ptr current = ClusterList<FS_TYPE>::read(params.current_filename, &cs);
    if (params.verbosity >= VerbosityNormal) stop_timer("done");

    // Perform tracking
    Tracking<FS_TYPE> tracking(params);
    tracking.track(previous, current);

    if (params.verbosity >= VerbosityNormal) {
        cout << "Resulting current list:" << endl;
        current->print();
    }

#if WITH_VTK
    if (params.write_vtk) {
        m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtu(current, cs, current->source_file);
        boost::filesystem::path path(params.current_filename);

        string modes_path = path.filename().stem().string() + "_modes.vtk";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_modes_vtk(modes_path, current->clusters, true);

        string centers_path = path.filename().stem().string() + "_centers.vtk";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_geometrical_cluster_centers_vtk(centers_path, current->clusters);
    }
#endif

    // Write results back
    if (params.verbosity >= VerbosityNormal) start_timer("-- Writing " + params.current_filename + " ... ");
    current->save_top_level_attributes();
    if (params.verbosity >= VerbosityNormal) stop_timer("done");

    // Clean up
    delete previous;
    delete current;
    return 0;
}
