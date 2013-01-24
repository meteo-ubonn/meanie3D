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
#include <netcdf>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;

/** Feature-space data type */
typedef double FS_TYPE;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm,
                        string &previous_filename,
                        string &current_filename,
                        string &tracking_variable_name,
                        Verbosity &verbosity)
{
    if ( vm.count("previous") == 0 )
    {
        cerr << "Missing previous cluster file argument --previous" << endl;
        
        exit( 1 );
    }
    
    try
    {
        previous_filename = vm["previous"].as<string>();
        
        NcFile previous_cluster_file( previous_filename, NcFile::read );
    }
    catch (const netCDF::exceptions::NcException &e )
    {
        cerr << "Exception opening file " << previous_filename << ":" << endl;
        
        cerr << e.what() << endl;
        
        exit(-1);
    }
    
    try
    {
        current_filename = vm["current"].as<string>();
        
        NcFile current_cluster_file( current_filename, NcFile::write );
    }
    catch (const netCDF::exceptions::NcException &e )
    {
        cerr << "Exception opening file " << current_filename << ":" << endl;
        
        cerr << e.what() << endl;
        
        exit(-1);
    }
    
    if ( vm.count("tracking-variable") > 0 )
    {
        tracking_variable_name = vm["tracking-variable"].as<string>();
    }
    
    // Verbosity
    
    unsigned short vb = vm["verbosity"].as<unsigned short>();
    
    if ( vb > VerbosityAll )
    {
        cerr << "Illegal value for parameter --verbosity. Only values from 0 .. 3 are allowed" << endl;
        
        exit( -1 );
    }
    else
    {
        verbosity = (Verbosity) vb;
    }
}

/**
 *
 *
 */
int main(int argc, char** argv)
{
    using namespace cfa::utils::timer;
    using namespace cfa::utils::visit;
    using namespace m3D;
    
    // Declare the supported options.
    
    program_options::options_description desc("Options");
    desc.add_options()
    ("help,h", "produce help message")
    ("previous,p", program_options::value<string>(), "Previous cluster file (netCDF)")
    ("current,c", program_options::value<string>(), "Current cluster file (netCDF)")
    ("tracking-variable,t", program_options::value<string>()->default_value("__default__"), "Variable used for histogram correlation. Defaults to the first variable that is not a dimension variable.")
    ("verbosity", program_options::value<unsigned short>()->default_value(1), "Verbosity level [0..3], 0=silent, 1=normal, 2=show details, 3=show all details). Default is 1.")
    ;
    
    program_options::variables_map vm;
    
    try
    {
        program_options::store( program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    }
    catch (std::exception &e)
    {
        cerr << "Error parsing command line: " << e.what() << endl;
        cerr << "Check meanie3D-track --help for command line options" << endl;
        exit(-1);
    }
    
    if ( vm.count("help")==1 || argc < 2 )
    {
        cout << desc << "\n";
        return 1;
    }
    
    // Evaluate user input
    
    string previous_filename, current_filename, tracking_variable_name;
    
    Verbosity verbosity;
    
    try
    {
        parse_commmandline(vm,previous_filename,current_filename,tracking_variable_name,verbosity);
    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    bool show_progress = (verbosity > VerbositySilent);
    
    // Read previous clusters
    
    ClusterList<FS_TYPE> prev_clusters;
    
    string prev_source, prev_parameters, prev_variable_names;
    
    vector<NcVar> previous_feature_variables, current_feature_variables;
    
    size_t prev_spatial_dims=0, curr_spatial_dims=0;
    
    NcFile *prev_file = NULL, *curr_file = NULL;
    
    ClusterList<FS_TYPE>::read(previous_filename,
                               &prev_file,
                               prev_clusters,
                               previous_feature_variables,
                               prev_spatial_dims,
                               prev_source,
                               prev_parameters,
                               prev_variable_names );
    
    // Read current clusters
    
    ClusterList<FS_TYPE> curr_clusters;
    
    string curr_source, curr_parameters, curr_variable_names;
    
    ClusterList<FS_TYPE>::read(current_filename,
                               &curr_file,
                               curr_clusters,
                               current_feature_variables,
                               curr_spatial_dims,
                               curr_source,
                               curr_parameters,
                               curr_variable_names );
    
    // Check if the feature variables match
    
    if ( previous_feature_variables != current_feature_variables )
    {
        cerr << "Incompatibe feature variables in the cluster files:" << endl;
        exit(-1);
    }
    
    // get the tracking variable
    
    NcVar tracking_var;
    
    if ( tracking_variable_name == "__default__" )
    {
        tracking_var = current_feature_variables[curr_spatial_dims];
    }
    else
    {
        bool found_tracking_var = false;
        
        for ( size_t i = 0; i < current_feature_variables.size(); i++ )
        {
            NcVar v = current_feature_variables[i];
            
            try {
                if ( v.getName() == tracking_variable_name )
                {
                    found_tracking_var = true;
                    
                    tracking_var = v;
                }
            }
            catch (const std::exception &e)
            {
                cerr << e.what() << endl;
                
                exit(-1);
            }
        }
        
        if ( !found_tracking_var)
        {
            cerr << "Tracking variable " << tracking_var.getName() << " is not part of the feature variables" << endl;

            exit(-1);
        }
    }
    
    // Perform tracking
    
    Tracking<FS_TYPE> tracking(1.0, 1.0, 1.0);

    tracking.track( prev_clusters, curr_clusters, current_feature_variables, tracking_var, curr_spatial_dims, verbosity );
    
    delete prev_file;
    
    delete curr_file;
    
    return 0;
};
