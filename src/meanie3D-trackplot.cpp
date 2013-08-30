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
#include <cf-algorithms/id.h>

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

typedef vector<size_t> bin_t;
typedef vector<string> svec_t;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm,
                        string &basename,
                        string &sourcepath,
                        svec_t &vtk_dim_names,
                        bin_t &length_histogram_bins,
                        bool &create_cluster_stats,
                        bin_t &cluster_histogram_bins,
                        bool &include_degenerates,
                        bool &create_cumulated_size_stats,
                        bin_t &size_histogram_bins,
                        bool &write_center_tracks_as_vtk,
                        bool &write_cumulated_tracks_as_vtk)
{
    if ( vm.count("basename") == 0 )
    {
        cerr << "Missing --basename argument" << endl;
        
        exit( 1 );
    }
    
    basename = vm["basename"].as<string>();
   
    sourcepath = vm["sourcepath"].as<string>();
    
    // Default is the dimensions
    
    if (vm.count("vtk-dimensions") > 0)
    {
        vtk_dim_names = vm["vtk-dimensions"].as<svec_t>();
    }
    
    include_degenerates = vm["exclude-degenerates"].as<bool>();
    
    length_histogram_bins = vm["length-histogram-classes"].as<bin_t>();

    create_cumulated_size_stats = vm["create-cumulated-size-statistics"].as<bool>();
    
    size_histogram_bins = vm["size-histogram-classes"].as<bin_t>();
    
    write_center_tracks_as_vtk = vm["write-center-tracks-as-vtk"].as<bool>();
    
    write_cumulated_tracks_as_vtk = vm["write-cumulated-tracks-as-vtk"].as<bool>();

    create_cluster_stats = vm["create-cluster-statistics"].as<bool>();
    
    cluster_histogram_bins = vm["cluster-histogram-classes"].as<bin_t>();
}

/** updates the size histogram
 */
void add_value_to_histogram(bin_t classes, bin_t &values, size_t size)
{
    bool found_bin = false;
    
    for (size_t i=0; i < classes.size() && !found_bin; i++)
    {
        if (i > 0)
        {
            if (size <= classes[i] && size > classes[i-1])
            {
                values[i] = values[i] + 1;
                found_bin = true;
            }
        }
        else
        {
            if (size <= classes[0])
            {
                values[i] = values[i] + 1;
                found_bin = true;
            }
        }
    }
    
    // not found? add to the last one
    
    if (!found_bin)
    {
        // Inform user
        
        cerr << "WARNING: value bigger than highest histogram class. Adding an 'over' bucket for values higher than " << values[classes.size()-1] << endl;
        
        // add the last class with value 0 if required
        
        if (values.size() == classes.size())
        {
            values.push_back(0);
        }
        
        // count up the 'over' bucket
        
        values[classes.size()] = values[classes.size()] + 1;
    }
}

/** Histogram output */

void print_histogram(bin_t classes, bin_t values, ofstream &file)
{
    for (size_t i=0; i < classes.size(); i++)
    {
        file << classes[i] << "," << values[i] << endl;
        cout << classes[i] << "," << values[i] << endl;
    }
    
    // 'Over' bucket?
    
    if (values.size() > classes.size())
    {
        file << " > " << classes[classes.size()-1] << "," << values.back() << endl;
        cout << " > " << classes[classes.size()-1] << "," << values.back() << endl;
    }
}

template <typename T>
double average(const vector<T> &v)
{
    // calculate the mean value
    double sum = 0.0;
    for (size_t i=0; i<v.size(); v++)
    {
        sum += boost::numeric_cast<double>(v[i]);
    }
    return sum / boost::numeric_cast<double>(v.size());
}

/** Takes the given 'distribution' and calculates the weighed average
 * @param distribution map
 * @param If <code>true</code> then the class '1' is ignored. 
 *        If <code>false</code> it counts normally.
 */
template <typename T>
double average(const map<size_t,size_t> &m, bool ignore_one = true)
{
    // calculate the mean value
    double sum = 0.0;
    
    size_t total_count = 0;
    
    map<size_t, size_t>::const_iterator mi;
    
    for (mi = m.begin(); mi != m.end(); mi++)
    {
        size_t key = mi->first;
        
        size_t val = mi->second;
        
        if (key == 1 && ignore_one)
        {
            continue;
        }
        
        total_count += val;
        
        sum += boost::numeric_cast<double>(mi->first * mi->second);
    }
    
    return sum / boost::numeric_cast<double>(total_count);
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
    
    size_t length_hist_default_values[] = {2,4,6,8,10,15,20,25,30,35,40,75,100,125,150,175,200};
    bin_t length_hist_default(length_hist_default_values,length_hist_default_values+17);
    
    size_t size_hist_default_values[] = {100,200,300,400,500,1000,2000,3000,4000,5000,7500,10000,15000,20000,25000,50000,100000};
    bin_t size_hist_default(size_hist_default_values,size_hist_default_values+17);
    
    size_t cluster_hist_default_values[] = {10,25,50,75,100,250,500,750,1000,1250,1500,1750,2000,2500,5000,7500,10000,15000,20000,50000,100000};
    bin_t cluster_hist_default(cluster_hist_default_values,cluster_hist_default_values+21);

    
    // Declare the supported options.
    
    program_options::options_description desc("Options");
    desc.add_options()
    ("help,h", "produce help message")
    ("basename,b", program_options::value<string>()->required(), "Previous cluster file (netCDF)")
    ("sourcepath,p", program_options::value<string>()->required()->default_value("."), "Current cluster file (netCDF)")
    ("length-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(length_hist_default),"List of track-length values for histogram bins")
    ("create-cluster-statistics",program_options::value<bool>()->default_value(true),"Evaluate each cluster in each track in terms of size.")
    ("cluster-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(cluster_hist_default),"List of cluster size values for histogram bins")
    ("exclude-degenerates", program_options::value<bool>()->default_value(true),"Exclude results of tracks of length one")
    ("create-cumulated-size-statistics",program_options::value<bool>()->default_value(true),"Evaluate each track in terms of cumulative size. Warning: this process takes a lot of memory.")
    ("size-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(size_hist_default),"List of cumulated track size values for histogram bins")
    ("write-center-tracks-as-vtk", program_options::value<bool>()->default_value(false),"Write tracks out as .vtk files")
    ("write-cumulated-tracks-as-vtk", program_options::value<bool>()->default_value(false),"Write cumulated tracks out as .vtk files. Only has effect if create-cumulated-size-statistics=true")
    ("vtk-dimensions", program_options::value<svec_t>()->multitoken(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
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
    bin_t length_histogram_classes;
    bin_t size_histogram_classes;
    bin_t cluster_histogram_classes;
    bool exclude_degenerates = true;
    bool write_center_tracks_as_vtk = false;
    bool write_cumulated_tracks_as_vtk = false;
    bool create_cumulated_size_statistics = false;
    bool create_cluster_statistics = false;
    svec_t vtk_dim_names;
    size_t number_of_degenerates = 0;
    
    try
    {
        parse_commmandline(vm,
                           basename,
                           sourcepath,
                           vtk_dim_names,
                           length_histogram_classes,
                           create_cluster_statistics,
                           cluster_histogram_classes,
                           exclude_degenerates,
                           create_cumulated_size_statistics,
                           size_histogram_classes,
                           write_center_tracks_as_vtk,
                           write_cumulated_tracks_as_vtk);
        
    }
    catch (const std::exception &e)
    {
        cerr << e.what() << endl;
        exit(-1);
    }

    // initialize histogram values
    
    bin_t length_histogram_values(length_histogram_classes.size(),0);
    bin_t size_histogram_values(size_histogram_classes.size(),0);
    bin_t cluster_histogram_values(cluster_histogram_classes.size(),0);
    
    // This map contains a mapping from the id type to lists of clusters.

    typename Tracking<T>::trackmap_t track_map;
    
    // Collect all files at sourcepath that start with basename and
    // end with -clusters.nc
    
    namespace fs = boost::filesystem;
    
    fs::path source_path(sourcepath);
    
    size_t spatial_dimensions = 0;
    
    std::string coords_filename;
    
    cout << "Collecting track data from NetCDF files: " << endl;
    
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
                    
                    cfa::id_t id = cluster->id;
                    
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
        
        cout << endl;
        cout << "Collating dimensions and dimension data ..." << endl;
        
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
        
        // Read "featurespace_dimensions"
        
        string fs_dimensions;
        coords_file->getAtt("featurespace_dimensions").getValues(fs_dimensions);
        vector<string> dim_names = from_string<string>(fs_dimensions);
        
        cout << "Attribute 'featurespace_dimensions' = " << dim_names << endl;
        
        ::m3D::utils::VisitUtils<T>::update_vtk_dimension_mapping(dim_names, vtk_dim_names);
        
        multimap<string,NcVar> vars = coords_file->getVars();
        multimap<string,NcVar>::iterator vi;
        
        // Only use dimensions, that have a variable of the
        // exact same name
        
        for (size_t di = 0; di < dim_names.size(); di++)
        {
            NcDim dim = coords_file->getDim(dim_names[di]);
            
            if (dim.isNull())
            {
                cerr << "ERROR: dimension " << dim_names[di] << " does not exist " << endl;
                exit(-1);
            }
            
            NcVar var = coords_file->getVar(dim_names[di]);
            
            if (var.isNull())
            {
                cerr << "ERROR: variable " << dim_names[di] << " does not exist " << endl;
                exit(-1);
            }
            
            dimensions.push_back(dim);
            dimension_vars.push_back(var);
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
        
        if (write_center_tracks_as_vtk)
        {
            ::m3D::utils::VisitUtils<T>::write_center_tracks_vtk(track_map, basename, spatial_dimensions, exclude_degenerates);
        }
        
        // Write cumulated tracks as netcdf and vtk files
        
        cout << "Cumulating track data: " << endl;
        
        // track length
        
        map<size_t,size_t> track_lengths;
        
        map<size_t,size_t> track_sizes;
        
        map<size_t,size_t> cluster_sizes;
        
        // Iterate over the collated tracks
        
        for (typename Tracking<T>::trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            size_t points_processed = 0;
            
            typename Tracking<T>::track_t *track = tmi->second;
            
            size_t track_length = track->size();
            
            if (track_length == 1)
            {
                number_of_degenerates++;
            }
            
            if (exclude_degenerates && track_length == 1)
            {
                continue;
            }
            
            cout << "Processing track #" << tmi->first << " (" << track->size() << " clusters)" << endl;
            
            // count up histogram
            
            add_value_to_histogram(length_histogram_classes, length_histogram_values, track_length);
            
            // add to distribution
            
            map<size_t,size_t>::iterator tlfi = track_lengths.find(track_length);
            
            if (tlfi==track_lengths.end())
            {
                track_lengths[track_length] = 1;
            }
            else
            {
                track_lengths[track_length] = (tlfi->second + 1);
            }
            
            // Skip the rest if cumulative size statistics and cluster
            // stats are both off
            
            if (!create_cumulated_size_statistics && !create_cluster_statistics)
            {
                continue;
            }

            // For each track, create an array index. Start with empty index.
            
            ArrayIndex<T> index(coord_system,false);
            
            // Iterate over the clusters in the track and sum up
            
            typename Tracking<T>::track_t::iterator ti;
            
            for (ti = track->begin(); ti != track->end(); ++ti)
            {
                Cluster<T> cluster = (*ti);
                
                if (create_cluster_statistics)
                {
                    size_t cluster_size = cluster.points.size();
                    
                    // add to histogram
                    
                    add_value_to_histogram(cluster_histogram_classes, cluster_histogram_values, cluster_size);
                    
                    // add to distribution
                    
                    map<size_t,size_t>::iterator csfi = cluster_sizes.find(cluster_size);
                    
                    if (csfi == cluster_sizes.end())
                    {
                        cluster_sizes[cluster_size] = 1;
                    }
                    else
                    {
                        cluster_sizes[cluster_size] = (csfi->second + 1);
                    }
                }
                
                // Skip the rest if cumulative size statistics are off
                
                if (!create_cumulated_size_statistics)
                {
                    continue;
                }
                
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
            
            // Update histogram
            
            add_value_to_histogram(size_histogram_classes, size_histogram_values, track_size);
            
            // Write out the cumulative file as netcdf and vtk
            
            // VTK
            
            if (write_cumulated_tracks_as_vtk)
            {
                string vtk_path = basename + "_cumulated_track_"+boost::lexical_cast<string>(tmi->first)+".vtk";
                
                ::cfa::meanshift::visit::VisitUtils<T>::write_pointlist_all_vars_vtk(vtk_path, &cumulatedList, vector<string>() );
            }
            
            // Don't need it anymore
            
            delete track;
            
            index.clear(true);
            
            // NETCDF
            
            cout <<"\t(processed " << points_processed << " points)" << endl;
            
        }
        
        cout << "done." << endl;
        
        // Write out statistical data in file file(s)
        
        string file_fn = basename+"_trackstats.txt";
        
        ofstream file(file_fn.c_str());
        
        // number of tracks
        
        file << "Tracking report for basename " + basename << endl;
        cout << "Tracking report for basename " + basename << endl;
        
        file << "(degenerates are excluded: " << (exclude_degenerates?"yes":"no") << ")" << endl;
        cout << "(degenerates are excluded: " << (exclude_degenerates?"yes":"no") << ")" << endl;
        
        file << endl;
        cout << endl;
        
        file << "Overall number of tracks: " << track_map.size() << endl;
        cout << "Overall number of tracks: " << track_map.size() << endl;
        
        file << "Number of degenerate tracks: " << number_of_degenerates << endl;
        cout << "Number of degenerate tracks: " << number_of_degenerates << endl;

        file << endl;
        cout << endl;
        
        //
        // length of tracks (in time steps)
        //
        
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;
        
        file << "Track length distribution" << endl;
        cout << "Track length distribution" << endl;
        
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<size_t>(track_lengths,true);
        file << "(Average track length = " << avg << ")" << endl;
        cout << "(Average track length = " << avg << ")" << endl;

        map<size_t, size_t>::reverse_iterator rend = track_lengths.rbegin();
        file << "(Maximum track length = " << rend->first << ")" << endl;
        cout << "(Maximum track length = " << rend->first << ")" << endl;
        
        file << "length,number" << endl;
        cout << "length,number" << endl;
        
        map<size_t,size_t>::iterator si;
        
        for (si=track_lengths.begin(); si!=track_lengths.end(); si++)
        {
            file << si->first << "," << si->second << endl;
            cout << si->first << "," << si->second << endl;        }
        
        file << endl;
        
        // Histogram
        
        file << endl;
        cout << endl;

        file << "Length Histogram:" << endl;
        cout << "Length Histogram:" << endl;

        file << "max length,number" << endl;
        cout << "max length,number" << endl;

        print_histogram(length_histogram_classes,length_histogram_values,file);
        
        
        //
        // Cluster sizes
        //
        
        file << endl;
        cout << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Cluster size distribution" << endl;
        cout << "Cluster size distribution" << endl;
        
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;
        
        avg = average<size_t>(cluster_sizes,true);
        file << "(Average cluster size = " << avg << ")" << endl;
        cout << "(Average cluster size = " << avg << ")" << endl;
        
        rend = cluster_sizes.rbegin();
        file << "(Maximum cluster size = " << rend->first << ")" << endl;
        cout << "(Maximum cluster size = " << rend->first << ")" << endl;
        
        file << "size,number" << endl;
        cout << "size,number" << endl;
        
        for (si=cluster_sizes.begin(); si!=cluster_sizes.end(); si++)
        {
            file << si->first << "," << si->second << endl;
            cout << si->first << "," << si->second << endl;
        }
        
        file << endl;
        cout << endl;
        
        // Histogram
        
        file << "Cluster Size Histogram:" << endl;
        cout << "Cluster Size Histogram:" << endl;
        
        file << "max size,number" << endl;
        cout << "max size,number" << endl;
        
        print_histogram(cluster_histogram_classes,cluster_histogram_values,file);


        //
        // Cumulative track stats
        //
        
        
        if (create_cumulated_size_statistics)
        {
            // size of tracks (in terms of cumulative pixels)
            
            file << endl;
            cout << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            file << "Track size distribution" << endl;
            cout << "Track size distribution" << endl;
            
            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            avg = average<size_t>(track_sizes,false);
            file << "(Average track size = " << avg << ")" << endl;
            cout << "(Average track size = " << avg << ")" << endl;
            
            rend = track_sizes.rbegin();
            file << "(Maximum track size = " << rend->first << ")" << endl;
            cout << "(Maximum track size = " << rend->first << ")" << endl;

            file << "size,number" << endl;
            cout << "size,number" << endl;
            
            for (si=track_sizes.begin(); si!=track_sizes.end(); si++)
            {
                file << si->first << "," << si->second << endl;
                cout << si->first << "," << si->second << endl;
            }
            
            file << endl;
            cout << endl;
            
            file << "Size Histogram:" << endl;
            cout << "Size Histogram:" << endl;
            
            file << "max size,number" << endl;
            cout << "max size,number" << endl;
            
            print_histogram(size_histogram_classes,size_histogram_values,file);

        }
        
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
