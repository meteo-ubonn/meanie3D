

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
                        map<int,double> &lower_thresholds,
                        map<int,double> &upper_thresholds,
                        double &scale,
                        int &weight_index,
                        string &parameters,
                        SearchParameters **search_params,
                        bool &write_vtk,
                        bool &write_weight_response,
                        vector<size_t> &vtk_dimension_indexes,
                        Verbosity &verbosity,
                        unsigned int &min_cluster_size,
                        PostAggregationMethod &post_aggregation,
                        double &drf_threshold,
                        vector<NcVar> &vtk_variables)
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
        NcVar var = file->getVar( *tok_iter );
        
        if (var.isNull())
        {
            cerr << "No variable '" << std::string(*tok_iter) << "' exists!" << endl;
            exit(-1);
        }
        
        variables.push_back( var );
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

    // Lower Thresholds
    
    if ( vm.count("lower-thresholds") > 0 )
    {
        boost::char_separator<char> equals("=");
        
        tokenizer tokens( vm["lower-thresholds"].as<string>(), sep );
        
        for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
        {
            std::string pair = *tok_iter;
            
            tokenizer subtokens( pair, equals );
            
            tokenizer::iterator subtoken_iter = subtokens.begin();
            
            std::string variableName = *subtoken_iter;
            
            NcVar variable;
            
            for (size_t i=0; i<variables.size(); i++)
            {
                if (variables[i].getName()==variableName)
                {
                    variable = variables[i];
                }
            }
            
            if (variable.isNull())
            {
                cerr << "No variable named " << variableName << " found. Check --lower-thresholds parameter" << endl;
                
                exit( 1 );
            }
            
            subtoken_iter++;
            
            if (subtoken_iter == subtokens.end())
            {
                cerr << "Missing threshold value for variable " << variableName << endl;
                
                exit( 1 );
            }
            
            const char* value = (*subtoken_iter).c_str();
            
            double doubleValue = strtod( value, (char **)NULL );
            
            lower_thresholds[variable.getId()] = doubleValue;
        }
    }
    
    // Upper Thresholds
    
    if ( vm.count("upper-thresholds") > 0 )
    {
        boost::char_separator<char> equals("=");
        
        tokenizer tokens( vm["upper-thresholds"].as<string>(), sep );
        
        for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
        {
            std::string pair = *tok_iter;
            
            tokenizer subtokens( pair, equals );
            
            tokenizer::iterator subtoken_iter = subtokens.begin();
            
            std::string variableName = *subtoken_iter;
            
            NcVar variable;
            
            for (size_t i=0; i<variables.size(); i++)
            {
                if (variables[i].getName()==variableName)
                {
                    variable = variables[i];
                }
            }
            
            if (variable.isNull())
            {
                cerr << "No variable named " << variableName << " found. Check --upper-thresholds parameter" << endl;
                
                exit( 1 );
            }
            
            subtoken_iter++;
            
            if (subtoken_iter == subtokens.end())
            {
                cerr << "Missing threshold value for variable " << variableName << endl;
                
                exit( 1 );
            }
            
            const char* value = (*subtoken_iter).c_str();
            
            double doubleValue = strtod( value, (char **)NULL );
            
            upper_thresholds[variable.getId()] = doubleValue;
        }
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
        weight_index = 0;
    }
    
    // VTK output?
    
    write_vtk = vm.count("write-clusters-as-vtk") > 0;
    
    write_weight_response = vm.count("write-cluster-weight-response") > 0;
    
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
    
    // vtk-variables
    
    if ( vm.count("write-variables-as-vtk") > 0 )
    {
        tokenizer bw_tokens( vm["write-variables-as-vtk"].as<string>(), sep );
        
        for ( tokenizer::iterator tok_iter = bw_tokens.begin(); tok_iter != bw_tokens.end(); ++tok_iter )
        {
            const char* bw = (*tok_iter).c_str();
            
            try
            {
                NcVar var = file->getVar(bw);
                
                if (var.isNull())
                {
                    cerr << "Can't open variable " << bw << " from NetCDF file. Check --write-variables-as-vtk" << endl;
                    exit(-1);
                }
                
                vtk_variables.push_back(var);
            }
            catch (const netCDF::exceptions::NcException &e)
            {
                cerr << "Can't find variable " << bw << " from NetCDF file. Check --write-variables-as-vtk" << endl;
                exit(-1);
            }
        }
    }
    
    if (vm.count("post-aggregate-clusters") > 0)
    {
        std::string val = vm["post-aggregate-clusters"].as<std::string>();
        
        if (val=="coalesce")
        {
            post_aggregation = PostAggregationMethodCoalescence;
        }
        else if (val == "drf")
        {
            post_aggregation = PostAggregationMethodDRF;
        }
        else if (val == "none")
        {
            post_aggregation = PostAggregationMethodNone;
        }
        else
        {
            cerr << "Unknown value for --post-aggregate-clusters (coalesce,drf,none)" << endl;
            exit(-1);
        }
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
    ("file,f", program_options::value<string>(), "CF-Metadata compliant NetCDF-file")
    ("output,o", program_options::value<string>(), "Name of output file for clustering results")
    ("dimensions,d", program_options::value<string>(), "Comma-separatred list of the dimensions to be used. The program expects dimension variables with identical names.")
    ("variables,v", program_options::value<string>(), "Comma-separated variables used to construct feature space. Do not include dimension variables")
    ("scale,s", program_options::value<double>()->default_value(NO_SCALE), "Scale parameter to pre-smooth the data with.")
    ("lower-thresholds", program_options::value<string>(), "Comma-separated list var1=val,var2=val,... of lower tresholds. Values below this are ignored when constructing feature space")
    ("upper-thresholds", program_options::value<string>(), "Comma-separated list var1=val,var2=val,... of lower tresholds. Values above this are ignored when constructing feature space")
    ("ranges,r", program_options::value<string>(), "Using this parameters means, you are choosing the range search method. Comma separated list of (dimensionless) ranges. Use in the order of (dim1,...dimN,var1,...,varN).")
    ("cluster-resolution,c",program_options::value<string>(), "When using KNN, specify the resolution of the clusters as comma separated list of (dimensionless) values. Use in the order of (dim1,...dimN,var1,...,varN).")
    ("post-aggregate-clusters,p",program_options::value<string>()->default_value("none"), "Run a post-processing step on raw meanshift clusters. Options are none (default),drf,coalesce")
    ("drf-threshold",program_options::value<double>()->default_value(0.65), "Used if -p=drf. DRF threshold for merging clusters. 1 means nothing is merged (=raw clusters). 0 means everything is merged as soon as it touches.")
    ("knn,k", program_options::value<long>(), "Using this parameter means, that you are choosing KNN-Method. The parameter value is the number of nearest neighbours to be used.")
    ("weight-variable,w", program_options::value<string>(), "variable to be used to weight meanshift. If none is given, the first variable is picked.")
    ("min-cluster-size,m",program_options::value<unsigned int>()->default_value(1u), "Keep only clusters of this minimum size at each pass")
    ("verbosity", program_options::value<unsigned short>()->default_value(1), "Verbosity level [0..3], 0=silent, 1=normal, 2=show details, 3=show all details). Default is 1.")
    ("write-clusters-as-vtk", "write clusters out in .vtk file format additionally (useful for visualization with visit for example)")
    ("write-cluster-weight-response","write out the clusters with weight responses as value")
    ("write-variables-as-vtk",program_options::value<string>(),"Comma separated list of variables that should be written out as VTK files (after applying scale/threshold)")
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
    map<int,double> lower_thresholds;   // ncvar.id / value
    map<int,double> upper_thresholds;   // ncvar.id / value
    vector<NcVar> vtk_variables;
    vector<double> cluster_resolution;
    int weight_index;
    string filename;
    string output_filename;
    string parameters;
    SearchParameters *search_params = NULL;
    bool write_vtk = false;
    bool write_weight_response = false;
    unsigned int min_cluster_size = 1;
    Verbosity verbosity = VerbosityNormal;
    double scale = NO_SCALE;
    PostAggregationMethod post_aggregation = PostAggregationMethodNone;
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
                           lower_thresholds,
                           upper_thresholds,
                           scale,
                           weight_index,
                           parameters,
                           &search_params,
                           write_vtk,
                           write_weight_response,
                           vtk_dimension_indexes,
                           verbosity,
                           min_cluster_size,
                           post_aggregation,
                           drf,
                           vtk_variables);
        
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
        
            cout << "\tcluster resolution =" << cluster_resolution << endl;
        }
        
        if ( min_cluster_size > 1 )
        {
            cout << "\tminimum cluster size = " << min_cluster_size << endl;
        }
        
        if ( scale != NO_SCALE )
        {
            double width = sqrt( ceil( -2.0*scale*log(0.01) ) ) / 2;
            
            cout << "\tpre-smoothing data with scale parameter " << scale << " ( kernel width = " << width << " )" << endl;
        }
        
        if ( !lower_thresholds.empty() )
        {
            cout << "\tusing lower thresholds " << vm["lower-thresholds"].as<string>() << endl;
        }

        if ( !upper_thresholds.empty() )
        {
            cout << "\tusing upper thresholds " << vm["upper-thresholds"].as<string>() << endl;
        }
        
        cout << "\tpost-processing of raw clusters: ";
        switch (post_aggregation)
        {
            case m3D::PostAggregationMethodNone:
                cout << "none.";
                break;
                
            case m3D::PostAggregationMethodCoalescence:
                cout << "coalescence";
                break;
                
            case m3D::PostAggregationMethodDRF:
                cout << "dynamic range factor = "<< drf;
                
            default:
                break;
        }
        
        cout << endl;

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
    //        PointIndex<FS_TYPE> *index = PointIndex<FS_TYPE>::create( fs );
    //        ClusterOperation<FS_TYPE> cop( fs, index );
    //        ClusterOperation<FS_TYPE>::reset_pass_counter();
    //        ClusterList<FS_TYPE> clusters = cop.cluster( search_params, cluster_resolution, kernel, weight_index, termcrit_epsilon, termcrit_iter, show_progress );
    //    }
    
    //
    // VISUALIZATION
    
    
    
    // Feature Space
    
    start_timer();
    
    // TODO: threshold should be applied as a filter somehow
    // this approach seems a little half-cocked
    
    FeatureSpace<FS_TYPE> *fs = new FeatureSpace<FS_TYPE>(filename,
                                                          coord_system,
                                                          variables,
                                                          lower_thresholds,
                                                          upper_thresholds,
                                                          show_progress );
#if WRITE_OFF_LIMITS_MASK
    fs->off_limits()->write("off_limits.vtk","off_limits");
#endif
    
    OASEWeightFunction<FS_TYPE> *weight_function = NULL;
    
    // Scale-Space smoothing
    
    if (scale != NO_SCALE)
    {
        // TODO: make decay a parameter or at least a constant
        
    	ScaleSpaceFilter<FS_TYPE> sf(scale, fs->coordinate_system->resolution(), 0.01, show_progress);
        
    	sf.apply(fs);
        
        if ( verbosity > VerbositySilent )
            cout << endl << "Constructing OASE weight function ...";
        
        weight_function = new OASEWeightFunction<FS_TYPE>(fs, sf.get_filtered_min(),sf.get_filtered_max());
        
        if ( verbosity > VerbositySilent )
            cout << " done." << endl;
        
        WeightThresholdFilter<FS_TYPE> wtf(weight_function,0.01, 10000, true);
        wtf.apply(fs);
    }
    else
    {
        if ( verbosity > VerbositySilent )
            cout << endl << "Constructing OASE weight function ...";
        
        weight_function = new OASEWeightFunction<FS_TYPE>(fs);
        
        if ( verbosity > VerbositySilent )
            cout << " done." << endl;
    }

    if (!vtk_variables.empty())
    {
        string filename_only = boost::filesystem::path(filename).filename().string();
        
        boost::filesystem::path destination_path = boost::filesystem::path(".");
        
        destination_path /= filename_only;
        
        //destination_path.replace_extension("vtk");
        
        destination_path.replace_extension();
        
        string dest_path = destination_path.generic_string();
        
        if ( verbosity > VerbositySilent )
            cout << "Writing featurespace-variables ...";
        
        cfa::utils::VisitUtils<FS_TYPE>::write_featurespace_variables_vtk(dest_path, fs, vtk_variables );
        
        if ( verbosity > VerbositySilent )
            cout << " done." << endl;
    }
    
    
#if WRITE_WEIGHT_FUNCTION
    boost::filesystem::path path(filename);
    std::string wfname = path.filename().stem().string() + "-weights.vtk";
    ::cfa::utils::VisitUtils<FS_TYPE>::write_weight_function_response(wfname, fs, weight_function);
#endif

    if ( verbosity == VerbosityAll )
        fs->print();
    
    // Calculate termcrit_epsilon if necessary
    
#if WRITE_MEANSHIFT_WEIGHTS
    vector<FS_TYPE> sample_point( 2 );
    sample_point[0] = -43.4621658325195;
    sample_point[1] = -4383.64453125;
    fs->weight_sample_points.push_back(sample_point);
#endif
    
    PointIndex<FS_TYPE> *index = PointIndex<FS_TYPE>::create( fs );
    
    //
    // Simple clustering
    //

    ClusterOperation<FS_TYPE> cop( fs, index );
    
    ClusterList<FS_TYPE> clusters = cop.cluster( search_params, kernel, weight_function, post_aggregation, drf, show_progress );
    
    // Sanity check
    
    // clusters.sanity_check( fs );
    
    if ( verbosity == VerbosityAll )
    {
        clusters.print();
    }
    
    // Axe weenies
    
    clusters.apply_size_threshold( min_cluster_size );
    
    // Number the result sequentially
    
    clusters.retag_identifiers();
    
    // Announce final results
    
    if ( verbosity > VerbositySilent )
    {
        cout << endl << "Final result: found " << clusters.clusters.size() << " objects: " << endl;
        
        clusters.print();
    }
    
    
    // Write out the cluster list
    
    if ( write_vtk && clusters.clusters.size() > 0)
    {
        boost::filesystem::path path(filename);
        ::m3D::utils::VisitUtils<FS_TYPE>::write_clusters_vtk( path.filename().string(), clusters.clusters, ranges, true );

        // MODES are needed for tagging with IDs
        
        string modes_path = path.filename().stem().string() + "-clusters_modes.vtk";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_modes_vtk( modes_path, clusters.clusters, true );
        
        string centers_path = path.filename().stem().string() + "-clusters_centers.vtk";
        ::m3D::utils::VisitUtils<FS_TYPE>::write_geometrical_cluster_centers_vtk( centers_path, clusters.clusters);
    }
    
    if ( write_weight_response && clusters.clusters.size() > 0 )
    {
        string wr_path = path.filename().stem().string() + "-clusters_weight";
        
        ::m3D::utils::VisitUtils<FS_TYPE>::write_cluster_weight_response_vtk(wr_path, clusters.clusters, weight_function, false);
    }
    
    if ( verbosity > VerbositySilent )
        cout << "Writing clusters to NetCDF file " << output_filename << " ..." << endl;
    
    clusters.write( output_filename );
    
    if ( verbosity > VerbositySilent )
        cout << "done." << endl;
    
    // mop up
    
    delete coord_system;
    
    delete kernel;
    
    delete index;
    
    delete fs;
    
    return 0;
};




