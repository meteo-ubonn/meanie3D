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

int main(int argc, char **argv) {
    // Declare the supported options.

    tracking_param_t params = Tracking<FS_TYPE>::defaultParams();

    program_options::options_description desc("Options");
    utils::add_standard_options(desc);
    add_tracking_options<FS_TYPE>(desc, params);

    program_options::variables_map vm;
    try {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    } catch (std::exception &e) {
        cerr << "FATAL:error parsing command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Get the command line content
    Verbosity verbosity;
    try {
        utils::get_standard_options(argc, vm, desc, verbosity);
        params.verbosity = verbosity;
        get_tracking_parameters<FS_TYPE>(vm, params);
    } catch (const std::exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (params.verbosity > VerbositySilent) {
        cout << "----------------------------------------------------" << endl;
        cout << "meanie3D-track" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << endl;
        print_tracking_params(params, vm);
    }

    // Read previous clusters
    if (params.verbosity >= VerbosityNormal) start_timer("Reading " + params.previous_filename + " ... ");
    ClusterList<FS_TYPE>::ptr previous = ClusterList<FS_TYPE>::read(params.previous_filename);
    if (params.verbosity >= VerbosityNormal) stop_timer("done");

#if WITH_VTK
    utils::set_vtk_dimensions_from_args<FS_TYPE>(vm, previous->dimensions);
#endif

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
    current->save();
    if (params.verbosity >= VerbosityNormal) stop_timer("done");

    // Clean up
    delete previous;
    delete current;

    return EXIT_SUCCESS;
}
