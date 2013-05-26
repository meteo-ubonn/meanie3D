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
typedef double T;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm,
                        string &basename,
                        string &sourcepath)
{
    if ( vm.count("basename") == 0 )
    {
        cerr << "Missing --basename argument" << endl;
        
        exit( 1 );
    }
    
    basename = vm["basename"].as<string>();
   
    sourcepath = vm["sourcepath"].as<string>();
}

/** reads the dimensions from command line parameter --vtk-dimensions
 * and updates the graphic output packages.
 */
void update_vtk_dimension_mapping(program_options::variables_map vm, const NcFile *file )
{
    // VTK dimension mapping
    
    if ( vm.count("vtk-dimensions") > 0 )
    {
        vector<size_t> vtk_dimension_indexes;
        
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
        
        ::cfa::utils::VisitUtils<T>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        
        ::m3D::utils::VisitUtils<T>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
        
        delete file;
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
    ("basename,b", program_options::value<string>(), "Previous cluster file (netCDF)")
    ("sourcepath,p", program_options::value<string>()->default_value("."), "Current cluster file (netCDF)")
    ("vtk-dimensions,d", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
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
        cerr << "Check meanie3D-trackplot --help for command line options" << endl;
        exit(-1);
    }
    
    if ( vm.count("help")==1 || argc < 2 )
    {
        cout << desc << "\n";
        return 1;
    }
    
    // Select the correct point factory
    PointFactory<T>::set_instance( new M3DPointFactory<T>() );
    
    // Evaluate user input
    
    string basename,sourcepath;
    
    try
    {
        parse_commmandline(vm, basename, sourcepath);
    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }
    
    // This map contains a mapping from the id type to lists of clusters.

    typename Tracking<T>::trackmap_t track_map;
    
    // Collect all files at sourcepath that start with basename and
    // end with -clusters.nc
    
    namespace fs = boost::filesystem;
    
    fs::path source_path(sourcepath);
    
    size_t spatial_dimensions = 0;
    
    std::string coords_filename;
    
    cout << "Collecing track data from NetCDF files: " << endl;
    
    if (fs::is_directory(source_path))
    {
        fs::directory_iterator dir_iter(source_path);
        fs::directory_iterator end;
        
        // Iterate over the "-clusters.nc" - files in sourcepath
        
        while (dir_iter != end)
        {
            fs::path f = dir_iter->path();
            
            if (fs::is_regular_file(f) && boost::algorithm::ends_with(f.filename().generic_string(), "-clusters.nc"))
            {
                // read the ClusterList from the file
                
                cout << "Processing " << f.filename().generic_string() << " ... ";
                
                // Remember the last one
                
                coords_filename = f.generic_string();
                
                typename ClusterList<T>::ptr cluster_list = ClusterList<T>::read( f.generic_string() );
                
                if (spatial_dimensions==0)
                {
                    spatial_dimensions = cluster_list->dimensions.size();
                }
                else if (spatial_dimensions != cluster_list->dimensions.size())
                {
                    cerr << "Spatial dimensions must remain identical across the track" << endl;
                    exit(-1);
                }
                
                // Iterate over the clusters in the list we just read
                
                typename Cluster<T>::list::iterator ci;
                
                for (ci=cluster_list->clusters.begin(); ci != cluster_list->clusters.end(); ++ci)
                {
                    typename Cluster<T>::ptr cluster = (*ci);
                    
                    cfa::meanshift::id_t id = cluster->id;
                    
                    typename Tracking<T>::trackmap_t::const_iterator ti = track_map.find(id);
                    
                    typename Tracking<T>::track_t * tm = NULL;
                    
                    if (ti==track_map.end())
                    {
                        // new entry
                        
                        tm = new Tracking<T>::track_t();
                        
                        track_map[id] = tm;
                    }
                    else
                    {
                        tm = ti->second;
                    }
                    
                    // append the cluster to the end of the track. Note that this
                    // is a copy operation. It's expensive, but it saves us from
                    // having to keep track of open files and stuff.
                    
                    tm->push_back( *cluster );
                }
                
                // Clean up
                
                delete cluster_list;
                
                cout << "done." << endl;
            }
            
            dir_iter++;
        }
        
        cout << endl;
        
        // Construct coordinate system
        
        CoordinateSystem<T> *coord_system = NULL;
        
        // This file should not be closed until the code has
        // run through, or netCDF will upchuck exceptions when
        // accessing coordinate system
        
        netCDF::NcFile *coords_file = NULL;

        // collate dimensions
        vector<NcDim> dimensions;
        vector<NcVar> dimension_vars;
        
        coords_file = new NcFile(coords_filename,NcFile::read);
        
        multimap<string,NcDim> dims = coords_file->getDims();
        multimap<string,NcDim>::iterator di;
        
        multimap<string,NcVar> vars = coords_file->getVars();
        multimap<string,NcVar>::iterator vi;
        
        // Only use dimensions, that have a variable of the
        // exact same name
        
        for (di=dims.begin(); di!=dims.end(); di++)
        {
            vi = vars.find(di->first);
            
            if (vi != vars.end())
            {
                dimensions.push_back(di->second);
                dimension_vars.push_back(vi->second);
            }
        }
        
        coord_system = new CoordinateSystem<T>(dimensions,dimension_vars);
        
        // Now we have a cluster map. Let's print it for debug purposes
        
        cout << "Keying up tracks: " << endl;
        
        for (typename Tracking<T>::trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            typename Tracking<T>::track_t *track = tmi->second;
            
            cout << "Track #" << tmi->first << " (" << track->size() << " clusters)" << endl;
            
            size_t i = 0;
            
            for ( typename Tracking<T>::track_t::iterator ti=track->begin(); ti!=track->end(); ti++ )
            {
                cout << "\t[" << i++ << "] x=" << ti->geometrical_center(spatial_dimensions) << endl;
            }
        }
        
        // Write center tracks
        
        ::m3D::utils::VisitUtils<T>::write_center_tracks_vtk(track_map, basename, spatial_dimensions);
        
        // Write cumulated tracks as netcdf and vtk files
        
        cout << "Cumulating track data: " << endl;
        
        // track length
        
        map<size_t,size_t> track_lengths;
        
        map<size_t,size_t> track_sizes;
        
        // Iterate over the collated tracks
        
        for (typename Tracking<T>::trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            size_t points_processed = 0;
            
            typename Tracking<T>::track_t *track = tmi->second;
            
            cout << "Processing track #" << tmi->first << " (" << track->size() << " clusters)" << endl;
    
            // For each track, create an array index. Start with empty index.
            
            ArrayIndex<T> index(coord_system);
            
            // Iterate over the clusters in the track and sum up
            
            typename Tracking<T>::track_t::iterator ti;
            
            for (ti = track->begin(); ti != track->end(); ++ti)
            {
                Cluster<T> cluster = (*ti);
                
                // Iterate over the points of the cluster
                
                typename Point<T>::list::iterator pi;
                
                for (pi = cluster.points.begin(); pi != cluster.points.end(); ++pi)
                {
                    Point<T>::ptr p = *pi;
                    
                    if (p->gridpoint.empty())
                    {
                        typename CoordinateSystem<T>::GridPoint gp = coord_system->newGridPoint();
                        
                        coord_system->reverse_lookup(p->coordinate, gp);
                        
                        p->gridpoint = gp;
                    }
                    
                    typename Point<T>::ptr indexed = index.get(p->gridpoint);
                    
                    if (indexed==NULL)
                    {
                        index.set(p->gridpoint, p);
                        
                        indexed = index.get(p->gridpoint);
                    }
                    
                    // only add up the value range
                    
                    for (size_t k = p->gridpoint.size(); k < indexed->values.size(); k++)
                    {
                        indexed->values[k] += p->values[k];
                    }
                    
                    points_processed++;
                }
            }
            
            // Extract the point list
            
            typename Point<T>::list cumulatedList;
            
            index.replace_points(cumulatedList);
            
            // Count the number of tracks with the given length (where
            // track length is the number of clusters it is comprised of)

            size_t track_length = track->size();
            
            map<size_t,size_t>::iterator tlfi = track_lengths.find(track_length);
            
            if (tlfi==track_lengths.end())
            {
                track_lengths[track_length] = 1;
            }
            else
            {
                track_lengths[track_length] = (tlfi->second + 1);
            }

            // Count the number of tracks with the given size, where
            // size is the number of points in it's cumulative list
            
            size_t track_size = cumulatedList.size();
            
            tlfi = track_sizes.find(track_size);
            
            if (tlfi==track_sizes.end())
            {
                track_sizes[track_size] = 1;
            }
            else
            {
                track_sizes[track_size] = (tlfi->second + 1);
            }
            
            // Write out the cumulative file as netcdf and vtk
            
            // VTK
            
            string vtk_path = basename + "_cumulated_track_"+boost::lexical_cast<string>(tmi->first)+".vtk";
            
            ::cfa::meanshift::visit::VisitUtils<T>::write_pointlist_all_vars_vtk(vtk_path, &cumulatedList, vector<string>() );
                               
            // Don't need it anymore
            
            delete track;
            
            // NETCDF
            
            cout <<"\t(processed " << points_processed << " points)" << endl;
            
        }
        
        cout << "done." << endl;
        
        // Write out statistical data in CSV file(s)
        
        string csv_fn = basename+"_trackstats.txt";
        
        ofstream csv(csv_fn.c_str());
        
        // number of tracks
        
        csv << "Tracking report for basename " + basename << endl;
        csv << endl;
        
        csv << "Overall number of tracks: " << track_map.size() << endl;
        csv << endl;
        
        // length of tracks (in time steps)
        
        csv << "length,number" << endl;
        
        map<size_t,size_t>::iterator si;
        
        for (si=track_lengths.begin(); si!=track_lengths.end(); si++)
        {
            csv << si->first << "," << si->second << endl;
        }
        
        csv << endl;
        
        // size of tracks (in terms of cumulative pixels)
        
        csv << "size,number" << endl;
        
        for (si=track_sizes.begin(); si!=track_sizes.end(); si++)
        {
            csv << si->first << "," << si->second << endl;
        }
        
        csv << endl;

        // close the cluster file we used to get the coordinate
        // system data
        
        delete coords_file;
    }
    else
    {
        cerr << "Argument --sourcepath does not point to a directory (sourcepath="+sourcepath+")" << endl;
        exit(-1);
    }
    
    
    return 0;
};
