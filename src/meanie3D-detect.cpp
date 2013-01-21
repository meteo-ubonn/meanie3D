

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

/** Feature-space data type */
typedef double FS_TYPE;

/** Verbosity */
typedef enum {
    VerbositySilent = 0,
    VerbosityNormal = 1,
    VerbosityDetails = 2,
    VerbosityAll = 3
} Verbosity;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline( program_options::variables_map vm,
                        NcFile **filePtr,
                        string &filename,
                        string &output_filename,
                        vector<NcDim> &dimensions,
                        vector<NcVar> &dimension_variables,
                        vector<NcVar> &variables,
                        map<NcVar,double> **thresholds,
                        double &scale,
                        int &weight_index,
                        string &parameters,
                        SearchParameters **search_params,
                        bool &write_vtk,
                        vector<size_t> &vtk_dimension_indexes,
                        Verbosity &verbosity,
                        unsigned int &min_cluster_size,
                        double &drf_threshold )
{
    if ( vm.count("file") == 0 )
    {
        cerr << "Missing input file argument" << endl;
        
        exit( 1 );
    }
    
    filename = vm["file"].as<string>();
    
    try
    {
        output_filename = vm["output"].as<string>();
    }
    catch ( const boost::exception& e )
    {
        cerr << "Missing parameter -o " << endl;
        
        exit(-1);
    }
    
    // Open NetCDF file
    
    NcFile *file = NULL;
    
    try
    {
        file = new NcFile( filename, NcFile::read );
    }
    catch (const netCDF::exceptions::NcException &e)
    {
        cerr << "ERROR opening file '" << filename << "' : " << e.what() << endl;
        exit(-1);
    }
    
    *filePtr = file;
    
    // Extract dimensions
    
    if ( vm.count("dimensions") == 0 )
    {
        cerr << "Missing parameter --dimensions" << endl;
        
        exit( 1 );
    }
    
    // parse dimension list
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    boost::char_separator<char> sep(",");
    
    tokenizer dim_tokens( vm["dimensions"].as<string>(), sep );
    
    for ( tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter )
    {
        const char* name = (*tok_iter).c_str();
        
        dimensions.push_back( file->getDim( name ) );
        
        dimension_variables.push_back( file->getVar( name ) );
    }
    
    parameters = parameters + "dimensions=" + vm["dimensions"].as<string>()+ " ";
    
    
    // parse variables
    
    if ( vm.count("variables") == 0 )
    {
        cerr << "Missing mandatory parameter --variables" << endl;
        
        exit( -1 );
    }
    
    tokenizer var_tokens( vm["variables"].as<string>(), sep );
    
    for ( tokenizer::iterator tok_iter = var_tokens.begin(); tok_iter != var_tokens.end(); ++tok_iter )
    {
        variables.push_back( file->getVar( *tok_iter ) );
    }
    
    parameters = parameters + "variables=" + vm["variables"].as<string>()+ " ";
    
    
    // Range-Search or KNN?
    
    if ( vm.count("knn") == 0 && vm.count("ranges") == 0 )
    {
        cerr << "Missing mandatory parameters --knn or --ranges" << endl;
        
        exit( -1 );
    }
    
    if ( vm.count("knn") > 0 && vm.count("ranges") > 0 )
    {
        cerr << "You must choose either --knn or --ranges" << endl;
        
        exit( -1 );
    }
    
    // parse ranges
    
    if ( vm.count("ranges") > 0 )
    {
        tokenizer bw_tokens( vm["ranges"].as<string>(), sep );
        
        vector<FS_TYPE> ranges;
        
        for ( tokenizer::iterator tok_iter = bw_tokens.begin(); tok_iter != bw_tokens.end(); ++tok_iter )
        {
            const char* bw = (*tok_iter).c_str();
            
            ranges.push_back( (FS_TYPE) strtod( bw, (char **)NULL ) );
        }
        
        if ( ranges.size() != dimension_variables.size() + variables.size() )
        {
            cerr << "Please provide " << dimension_variables.size() + variables.size() << " bandwidth values" << endl;
            
            exit( 1 );
        }
        
        parameters = parameters + "ranges=" + vm["ranges"].as<string>();
        
        RangeSearchParams<FS_TYPE> *p = new RangeSearchParams<FS_TYPE>( ranges );
        
        *search_params = p;
    }
    else if ( vm.count("knn") > 0 )
    {
        // cluster_resolution
        
        if ( vm.count("cluster-resolution") == 0 )
        {
            cerr << "When using knn, you must specify --cluster-resolution " << endl;
            exit( -1 );
        }
        
        size_t k = vm["knn"].as<long>();
        
        vector<FS_TYPE> cluster_resolution;
        
        tokenizer tokens( vm["cluster-resolution"].as<string>(), sep );
        
        for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
        {
            const char* value = (*tok_iter).c_str();
            
            cluster_resolution.push_back( strtod( value, (char **)NULL ) );
        }
        
        if ( cluster_resolution.size() != ( dimension_variables.size() + variables.size() ) )
        {
            cerr << "Please provide " << dimension_variables.size() + variables.size() << " cluster-resolution values" << endl;
            
            exit( 1 );
        }
        
        KNNSearchParams<FS_TYPE> *p = new KNNSearchParams<FS_TYPE>( k,cluster_resolution );
        
        *search_params = p;
        
        parameters = parameters + "knn=" + boost::lexical_cast<string>(k) + ",resolution=" + vm["cluster-resolution"].as<string>();
    }
    
    if ( vm.count("thresholds") > 0 )
    {
        *thresholds = FeatureSpace<FS_TYPE>::NO_THRESHOLDS;
        
        map<NcVar,double> *th = new map<NcVar,double>();
        
        tokenizer tokens( vm["thresholds"].as<string>(), sep );
        
        size_t var_index = 0;
        
        for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
        {
            if ( var_index > (variables.size()-1) )
            {
                cerr << "Please provide " << variables.size() << " thresholds values" << endl;
                
                exit( 1 );
            }
            
            const char* value = (*tok_iter).c_str();
            
            NcVar var = variables[var_index];
            
            th->operator[](var) = strtod( value, (char **)NULL );
            
            var_index++;
        }
        
        if ( th->size() != variables.size() )
        {
            cerr << "Please provide " << variables.size() << " thresholds values" << endl;
            
            exit( 1 );
        }
        
        *thresholds = th;
    }
    
    scale = vm["scale"].as<double>();
    
    drf_threshold = vm["drf-threshold"].as<double>();
    
    // get weighing variable
    
    if ( vm.count("weight-variable") > 0 )
    {
        string var_name = vm["weight-variable"].as<string>();
        
        // test if the variable exists
        
        NcVar var = file->getVar( var_name );
        
        if ( var.isNull() )
        {
            cerr << "No such variable " << var_name << endl;
            
            exit( -1 );
        }
        
        // test if the variable is in the list of variables
        
        int var_index = -1;
        
        for ( size_t index=0; index < variables.size(); index++ )
        {
            if ( variables[index] == var )
            {
                var_index = index;
                
                break;
            }
        }
        
        if ( var_index < 0 )
        {
            cerr << "Weight variable must be in the list of variables" << endl;
            
            exit(-1);
        }
        
        // figure out the weight_index
        
        weight_index = var_index;
        
        parameters = parameters + " weight-variable=" + var_name;
    }
    else
    {
        weight_index = MeanshiftOperation<FS_TYPE>::NO_WEIGHT;
    }
    
    // VTK output?
    
    write_vtk = vm.count("write-clusters-as-vtk") > 0;
    
    // VTK dimension mapping
    
    if ( vm.count("vtk-dimensions") > 0 )
    {
        // parse dimension list
        
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        
        boost::char_separator<char> sep(",");
        
        string str_value = vm["vtk-dimensions"].as<string>();
        
        tokenizer dim_tokens( str_value, sep );
        
        for ( tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter )
        {
            const char* name = (*tok_iter).c_str();
            
            NcDim dim = file->getDim( name );
            
            vector<NcDim>::const_iterator fi = find( dimensions.begin(), dimensions.end(), dim );
            
            if ( fi == dimensions.end() )
            {
                cerr << "--vtk-dimension parameter " << dim.getName() << " is not part of --dimensions" << endl;
                exit(-1);
            }
            
            size_t index = fi - dimensions.begin();
            
            vtk_dimension_indexes.push_back( index );
        }
        
        if ( vtk_dimension_indexes.size() != dimensions.size() )
        {
            cerr << "The number of vtk-dimensions must be identical to dimensions" << endl;
            
            exit(-1);
        }
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
    
    min_cluster_size = vm["min-cluster-size"].as<unsigned int>();
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
    ("file,f", program_options::value<string>(), "CF-Metadata compliant NetCDF-file")
    ("output,o", program_options::value<string>(), "Name of output file for clustering results")
    ("dimensions,d", program_options::value<string>(), "Comma-separatred list of the dimensions to be used. The program expects dimension variables with identical names.")
    ("variables,v", program_options::value<string>(), "Comma-separated variables used to construct feature space. Do not include dimension variables")
    ("scale,s", program_options::value<double>()->default_value(NO_SCALE), "Scale parameter to pre-smooth the data with.")
    ("thresholds,t", program_options::value<string>(), "Comma-separated list of tresholds, one per variable. Values below this are ignored when constructing feature space")
    ("ranges,r", program_options::value<string>(), "Using this parameters means, you are choosing the range search method. Comma separated list of (dimensionless) ranges. Use in the order of (dim1,...dimN,var1,...,varN).")
    ("cluster-resolution,c",program_options::value<string>(), "When using KNN, specify the resolution of the clusters as comma separated list of (dimensionless) values. Use in the order of (dim1,...dimN,var1,...,varN).")
    ("drf-threshold",program_options::value<double>()->default_value(0.65), "Dynamic range factor threshold for merging clusters. 1 means nothing is merged (=raw clusters). 0 means everything is merged as soon as it touches.")
    ("knn,k", program_options::value<long>(), "Using this parameter means, that you are choosing KNN-Method. The parameter value is the number of nearest neighbours to be used.")
    ("weight-variable,w", program_options::value<string>(), "variable to be used to weigh meanshift")
    ("min-cluster-size,m",program_options::value<unsigned int>()->default_value(1u), "Keep only clusters of this minimum size at each pass")
    ("verbosity", program_options::value<unsigned short>()->default_value(1), "Verbosity level [0..3], 0=silent, 1=normal, 2=show details, 3=show all details). Default is 1.")
    ("write-clusters-as-vtk", "write clusters out in .vtk file format additionally (useful for visualization with visit for example)")
    ("vtk-dimensions", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
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
        cerr << "Check meanie3D-detect --help for command line options" << endl;
        exit(-1);
    }
    
    if ( vm.count("help")==1 || argc < 2 )
    {
        cout << desc << "\n";
        return 1;
    }
    
    // Evaluate user input
    
    NcFile *file = NULL;
    vector<NcDim> dimensions;
    vector<size_t> vtk_dimension_indexes;
    vector<NcVar> dimension_variables;
    vector<NcVar> variables;
    vector<double> ranges;
    map<NcVar,double> *thresholds = FeatureSpace<FS_TYPE>::NO_THRESHOLDS;
    vector<double> cluster_resolution;
    int weight_index;
    string filename;
    string output_filename;
    string parameters;
    SearchParameters *search_params = NULL;
    bool write_vtk = false;
    unsigned int min_cluster_size = 1;
    Verbosity verbosity = VerbosityNormal;
    double scale = NO_SCALE;
    double drf = 0.95;
    
    try
    {
        parse_commmandline( vm,
                           &file,
                           filename,
                           output_filename,
                           dimensions,
                           dimension_variables,
                           variables,
                           &thresholds,
                           scale,
                           weight_index,
                           parameters,
                           &search_params,
                           write_vtk,
                           vtk_dimension_indexes,
                           verbosity,
                           min_cluster_size,
                           drf);
        
        // Make the mapping known to the visualization routines
        
        ::cfa::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        ::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        
        // Select the correct point factory
        PointFactory<FS_TYPE>::set_instance( new M3DPointFactory<FS_TYPE>() );
        
    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    bool show_progress = (verbosity > VerbositySilent);
    
    
    // DEBUG
    //
    
    //    vector<double> scales;
    //    scales.push_back(0);
    //    scales.push_back(1);
    //    scales.push_back(10);
    //    scales.push_back(100);
    //    scales.push_back(1000);
    //    scales.push_back(5000);
    //
    //    CoordinateSystem<FS_TYPE> *cs = new CoordinateSystem<FS_TYPE>( dimensions, dimension_variables );
    //
    //    vector<double>::const_iterator si;
    //    for ( si=scales.begin(); si!=scales.end(); si++ )
    //    {
    //        double t = *si;
    //
    //        // Feature Space
    //        FeatureSpace<double> *sfs = new FeatureSpace<double>( filename, cs, variables, cluster_resolution, thresholds, t, show_progress );
    //
    //        delete sfs;
    //    }
    //
    //    delete cs;
    
    //
    // DEBUG
    
    //    // DEBUG
    //    //
    //
    //#include <cf-algorithms/utils/console_utils.h>
    //
    //    Meanshift::utils::ConsoleSpinner spinner;
    //
    //    spinner.start();
    //
    //    sleep( 10000 );
    //
    //    spinner.stop();
    //
    //    //
    //    // DEBUG
    
    
    
    RangeSearchParams<FS_TYPE> *p = dynamic_cast<RangeSearchParams<FS_TYPE> *>(search_params);
    
    if ( verbosity > VerbositySilent )
    {
        cout << "clustering is going to run with the following parameters:" << endl;
        
        cout << "\tinput file = " << filename << endl;
        
        cout << "\tdimensions = " << vm["dimensions"].as<string>() << endl;
        
        cout << "\tvariables = " << vm["variables"].as<string>() << endl;
        
        if ( p != NULL )
        {
            cout << "\tranges = " << p->bandwidth << endl;
        }
        else
        {
            KNNSearchParams<FS_TYPE> *knn = dynamic_cast<KNNSearchParams<FS_TYPE> *>(search_params);
            
            cout << "\tknn = " << knn->k << ", resolution = " << knn->resolution << endl;
        }
        
        cout << "\tcluster resolution =" << cluster_resolution << endl;
        
        if ( min_cluster_size > 1 )
        {
            cout << "\tminimum cluster size = " << min_cluster_size << endl;
        }
        
        if ( scale != NO_SCALE )
        {
            double width = sqrt( ceil( -2.0*scale*log(0.01) ) ) / 2;
            
            cout << "\tpre-smoothing data with scale parameter " << scale << " ( kernel width = " << width << " )" << endl;
        }
        
        if ( thresholds != FeatureSpace<FS_TYPE>::NO_THRESHOLDS )
        {
            cout << "\tusing thresholds " << vm["thresholds"].as<string>() << endl;
        }
        
        cout << "\tdynamic range factor = " << drf << endl;
        
        cout << "\toutput written to file: " << output_filename << endl;
        
        cout << "\tclusters written as vtk: " << (write_vtk ? "yes":"no") << endl;
    }
    
    // Construct Featurespace
    
    // Coordinate system
    
    CoordinateSystem<FS_TYPE> *coord_system = new CoordinateSystem<FS_TYPE>( dimensions, dimension_variables );
    Kernel<FS_TYPE> *kernel = new UniformKernel<FS_TYPE>(1.0);
    
    // VISUALIZAION
    //
    
    //    vector<double> scales;
    //    scales.push_back(0);
    //    scales.push_back(2);
    //    scales.push_back(4);
    //    scales.push_back(8);
    //    scales.push_back(16);
    //    scales.push_back(32);
    //    scales.push_back(64);
    //    scales.push_back(128);
    //    scales.push_back(256);
    //    scales.push_back(512);
    //
    //    cout << "Running analysis on scales " << scales << endl;
    //
    //    for ( size_t i=0; i<scales.size(); i++)
    //    {
    //        cout << endl << "== SCALE " << scales[i] << " ==" << endl;
    //        FeatureSpace<FS_TYPE> *fs = new FeatureSpace<FS_TYPE>( filename, coord_system, variables, cluster_resolution, thresholds, scales[i], show_progress );
    //        FeatureSpaceIndex<FS_TYPE> *index = FeatureSpaceIndex<FS_TYPE>::create( fs );
    //        ClusterOperation<FS_TYPE> cop( fs, index );
    //        ClusterOperation<FS_TYPE>::reset_pass_counter();
    //        ClusterList<FS_TYPE> clusters = cop.cluster( search_params, cluster_resolution, kernel, weight_index, termcrit_epsilon, termcrit_iter, show_progress );
    //    }
    
    //
    // VISUALIZATION
    
    
    
    // Feature Space
    
    start_timer();
    
    FeatureSpace<FS_TYPE> *fs = new FeatureSpace<FS_TYPE>( filename, coord_system, variables, thresholds, show_progress );
    
    // Scale-Space smoothing
    
    if (scale != NO_SCALE) {
    	ScaleSpaceFilter<FS_TYPE> sf(scale,show_progress);
        
    	sf.apply(fs);
    }
    
    if ( verbosity == VerbosityAll )
        fs->print();
    
    // Calculate termcrit_epsilon if necessary
    
#if WRITE_MEANSHIFT_WEIGHTS
    vector<FS_TYPE> sample_point( 2 );
    sample_point[0] = -43.4621658325195;
    sample_point[1] = -4383.64453125;
    fs->weight_sample_points.push_back(sample_point);
#endif
    
    FeatureSpaceIndex<FS_TYPE> *index = FeatureSpaceIndex<FS_TYPE>::create( fs );
    
    //
    // Simple clustering
    //
    
    ClusterOperation<FS_TYPE> cop( fs, index );
    
    ClusterList<FS_TYPE> clusters = cop.cluster( search_params, kernel, weight_index, drf, show_progress );
    
    // Sanity check
    
    clusters.sanity_check( fs );
    
    if ( verbosity == VerbosityAll )
    {
        clusters.print();
    }
    
    // Axe weenies
    
    clusters.apply_size_threshold( min_cluster_size );
    
    // Announce final results
    
    if ( verbosity > VerbositySilent )
    {
        cout << endl << "Final result: found " << clusters.clusters.size() << " objects: " << endl;
        
        clusters.print();
    }
    
    
    // Write out the cluster list
    
    if ( write_vtk )
    {
        boost::filesystem::path path(filename);
        
        ::m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtk( path.filename().string(), clusters.clusters, ranges, true );
    }
    
    if ( verbosity > VerbositySilent )
        cout << "Writing clusters to NetCDF file " << output_filename << " ..." << endl;
    
    clusters.write( output_filename, fs, parameters );
    
    if ( verbosity > VerbositySilent )
        cout << "done." << endl;
    
    // mop up
    
    delete coord_system;
    
    delete kernel;
    
    delete index;
    
    delete fs;
    
    return 0;
};




