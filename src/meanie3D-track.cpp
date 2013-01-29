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
                        bool &write_vtk,
                        vector<size_t> &vtk_dimension_indexes,
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
    
    // VTK dimension mapping
    
    if ( vm.count("vtk-dimensions") > 0 )
    {
        
        // Open the file for reading once more
        
        NcFile *file = new NcFile( current_filename, NcFile::read );
        
        
        // Read "featurespace_dimensions"
        
        string fs_dimensions;

        file->getAtt("featurespace_dimensions").getValues(fs_dimensions);

        vector<string> fs_dim_names = from_string<string>(fs_dimensions);
        
        
        // parse dimension list
        
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        
        boost::char_separator<char> sep(",");
        
        string str_value = vm["vtk-dimensions"].as<string>();
        
        tokenizer dim_tokens( str_value, sep );

        for ( tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter )
        {
            string name = *tok_iter;
            
            int index = index_of_first( fs_dim_names, name );
            
            if ( index < 0 )
            {
                cerr << "Invalid dimension '" << name << "'. Check parameter --vtk-dimensions" << endl;
                
                exit(-1);
            }
            
            vtk_dimension_indexes.push_back( (size_t)index );
        }
        
        delete file;
    }
    
    // --write-vtk
    
    write_vtk = vm.count("write-vtk") > 0;
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
    ("write-vtk,k","Write out the clusters as .vtk files for visit")
    ("vtk-dimensions", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
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
    
    // Select the correct point factory
    PointFactory<FS_TYPE>::set_instance( new M3DPointFactory<FS_TYPE>() );
    
    // Evaluate user input
    
    string previous_filename, current_filename, tracking_variable_name;
    Verbosity verbosity;
    bool write_vtk = false;
    vector<size_t> vtk_dimension_indexes;
    
    try
    {
        parse_commmandline(vm,previous_filename,current_filename,tracking_variable_name,write_vtk,vtk_dimension_indexes,verbosity);
        
        ::cfa::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        ::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;

    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    bool show_progress = (verbosity > VerbositySilent);
    
    // Read previous clusters
    
    typename ClusterList<FS_TYPE>::ptr previous = ClusterList<FS_TYPE>::read( previous_filename );

    // Read current clusters
    
    typename ClusterList<FS_TYPE>::ptr current = ClusterList<FS_TYPE>::read( current_filename );
    
    // Check if the feature variables match
    
    if ( previous->feature_variables != current->feature_variables )
    {
        cerr << "Incompatibe feature variables in the cluster files:" << endl;
        exit(-1);
    }
    
    // get the tracking variable
    
    NcVar tracking_var;
    
    if ( tracking_variable_name == "__default__" )
    {
        tracking_var = current->feature_variables[current->dimensions.size()];
    }
    else
    {
        bool found_tracking_var = false;
        
        for ( size_t i = 0; i < current->feature_variables.size(); i++ )
        {
            NcVar v = current->feature_variables[i];
            
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
    
    Tracking<FS_TYPE> tracking;
    
    tracking.track( previous, current, tracking_var, verbosity );
    
    // Write results back
    
    current->write( current_filename );
    
    if ( write_vtk )
    {
        m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtk( previous );
        
        m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtk( current );
    }
    
    // Clean up
    
    delete previous;
    
    delete current;
    
    return 0;
};