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
#include <boost/algorithm/string.hpp>

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
#include <unistd.h>

#include <meanie3D/meanie3D.h>
#include <radolan/radolan.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;
using namespace Radolan;

#pragma mark -
#pragma mark type definitions

/** Feature-space data type */
typedef double FS_TYPE;

typedef vector<float> fvec_t;
typedef vector<size_t> bin_t;
typedef vector<string> svec_t;

#pragma mark -
#pragma mark comand line parsing

void parse_commmandline(program_options::variables_map vm,
        string &basename,
        string &sourcepath,
        svec_t &vtk_dim_names,
        bool &create_length_stats,
        bin_t &length_histogram_bins,
        bool &create_speed_stats,
        fvec_t &speed_histogram_bins,
        bool & create_direction_stats,
        fvec_t &direction_histogram_classes,
        bool &create_cluster_stats,
        bin_t &cluster_histogram_bins,
        bool &include_degenerates,
        bool &create_cumulated_size_stats,
        bin_t &size_histogram_bins,
        bool &write_center_tracks_as_vtk,
        bool &write_cumulated_tracks_as_vtk,
        bool &write_gnuplot_files,
        bool &write_track_dictionary)
{
    // Version

    sourcepath = vm["sourcepath"].as<string>();

    if (vm.count("basename") == 0)
    {
        basename = ::m3D::utils::common_component(sourcepath, ".nc");
    } else
    {
        basename = vm["basename"].as<string>();
    }

    write_center_tracks_as_vtk = vm.count("write-center-tracks-as-vtk") > 0;

    write_gnuplot_files = vm.count("write-gnuplot-files") > 0;

    write_cumulated_tracks_as_vtk = vm.count("write-cumulated-tracks-as-vtk") > 0;

    write_track_dictionary = vm.count("write-track-dictionary") > 0;

    include_degenerates = vm["exclude-degenerates"].as<bool>();

    // Default is the dimensions

    if (vm.count("vtk-dimensions") > 0) {
        std::string s = vm["vtk-dimensions"].as<std::string>();
        if (contains(s,",") || contains(s,"(") || contains(s,")")) {
            cerr << "wrong arguments for --vtk-dimension. Example: --vtk-dimensions x y z" << endl;
            exit(EXIT_FAILURE);
        }
        vtk_dim_names = vm["vtk-dimensions"].as<svec_t>();
    }

    create_length_stats = vm.count("create-length-statistics") > 0;
    length_histogram_bins = vm["length-histogram-classes"].as<bin_t>();

    create_speed_stats = vm.count("create-speed-statistics") > 0;
    speed_histogram_bins = vm["speed-histogram-classes"].as<fvec_t>();

    create_direction_stats = vm.count("create-direction-statistics") > 0;
    direction_histogram_classes = vm["direction-histogram-classes"].as<fvec_t>();

    create_cumulated_size_stats = vm.count("create-cumulated-size-statistics") > 0;
    size_histogram_bins = vm["size-histogram-classes"].as<bin_t>();

    create_cluster_stats = vm.count("create-cluster-statistics") > 0;
    cluster_histogram_bins = vm["cluster-histogram-classes"].as<bin_t>();
}

#pragma mark -
#pragma mark histogram helpers

/** updates the size histogram
 * @param classes
 * @param counts
 * @param value
 * @param exceeded_max_class contains <code>true</code> if the
 *        value was added to the 'over' bucket. <code>false</code> else
 */
template <typename T>
void add_value_to_histogram(const vector<T> &classes,
        bin_t &counts,
        T value,
        bool &exceeded_max_class)
{
    bool found_bin = false;

    exceeded_max_class = false;

    for (size_t i = 0; i < classes.size() && !found_bin; i++)
    {
        if (i > 0)
        {
            if (value <= classes[i] && value > classes[i - 1])
            {
                counts[i] = counts[i] + 1;
                found_bin = true;
            }
        } else
        {
            if (value <= classes[0])
            {
                counts[i] = counts[i] + 1;
                found_bin = true;
            }
        }
    }

    // not found? add to the last one

    if (!found_bin)
    {
        // add the last class with value 0 if required

        if (counts.size() == classes.size())
        {
            counts.push_back(0);
            cerr << "WARNING: value bigger than highest histogram class. Adding an 'over' bucket for values higher than " << classes[classes.size() - 1] << endl;
        }

        // count up the 'over' bucket

        counts[classes.size()] = counts[classes.size()] + 1;

        exceeded_max_class = true;
    }
}

/** Histogram output */

template <typename T>
void print_histogram(const vector<T> &classes,
        const bin_t &values,
        ofstream &file)
{
    for (size_t i = 0; i < classes.size(); i++)
    {
        file << classes[i] << "," << values[i] << endl;
        cout << classes[i] << "," << values[i] << endl;
    }

    // 'Over' bucket?

    if (values.size() > classes.size())
    {
        file << " > " << classes[classes.size() - 1] << "," << values.back() << endl;
        cout << " > " << classes[classes.size() - 1] << "," << values.back() << endl;
    }
}

template <typename T>
void write_histogram(const std::string &filename,
        const std::string &x,
        const std::string &y,
        const vector<T> &classes,
        const bin_t &values)
{
    ofstream file(filename.c_str());

    file << "#" << x << "  " << y << endl;

    for (size_t i = 0; i < classes.size(); i++)
    {
        file << classes[i] << "  " << values[i] << endl;
    }

    file.close();
}

template <typename T>
void write_values(const std::string &filename,
        const std::string &x,
        const std::string &y,
        const map<T, size_t> &values)
{
    ofstream file(filename.c_str());

    file << "track length  number of tracks" << endl;

    typename map<T, size_t>::const_iterator si;

    for (si = values.begin(); si != values.end(); si++)
    {
        file << si->first << "  " << si->second << endl;
    }

    file.close();
}

template <typename T>
double average(const vector<FS_TYPE> &v)
{
    // calculate the mean value
    double sum = 0.0;
    for (size_t i = 0; i < v.size(); i++)
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
double average(const map<T, size_t> &m, bool ignore_one = true)
{
    // calculate the mean value
    double sum = 0.0;

    size_t total_count = 0;

    typename std::map<T, size_t>::const_iterator mi;

    for (mi = m.begin(); mi != m.end(); mi++)
    {
        T key = mi->first;

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

#pragma mark -
#pragma mark Memory management

#if APPLE 
#include <mach/vm_statistics.h>
#include <mach/mach_types.h> 
#include <mach/mach_init.h>
#include <mach/mach_host.h>

/** Calculates available and used memory in megabytes
 * 
 * @param free_memory
 * @param used_memory
 */
void get_memory(double &free_memory, double &used_memory)
{
    const double MB = 1024.0 * 1024.0;

    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics_data_t vm_stats;

    mach_port = mach_host_self();

    count = sizeof (vm_stats) / sizeof (natural_t);

    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
            KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
            (host_info_t) & vm_stats, &count))
    {
        free_memory = ((double) vm_stats.free_count * (int64_t) page_size) / MB;



        used_memory = (((double) vm_stats.active_count +
                (double) vm_stats.inactive_count +
                (double) vm_stats.wire_count) * (double) page_size) / MB;
    }
}
#endif

#pragma mark -
#pragma mark main

/**
 *
 *
 */
int main(int argc, char** argv)
{
    using namespace m3D;
    using namespace m3D::utils::vectors;

    size_t length_hist_default_values[] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110};
    bin_t length_hist_default(length_hist_default_values, length_hist_default_values + 22);

    size_t size_hist_default_values[] = {100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000, 20000, 25000, 50000, 100000};
    bin_t size_hist_default(size_hist_default_values, size_hist_default_values + 17);

    size_t cluster_hist_default_values[] = {10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 10000000};
    bin_t cluster_hist_default(cluster_hist_default_values, cluster_hist_default_values + 11);

    float speed_hist_default_values[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60};
    fvec_t speed_hist_default(speed_hist_default_values, speed_hist_default_values + 30);

    float direction_hist_default_values[] = {15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 175, 200, 215, 230, 245, 260, 275, 300, 315, 330, 345, 360};
    fvec_t direction_hist_default(direction_hist_default_values, direction_hist_default_values + 22);

    // Declare the supported options.

    program_options::options_description desc("Options");
    desc.add_options()
            ("help,h", "produce help message")
            ("version", "print version information and exit")
            ("basename,b", program_options::value<string>()->implicit_value(""), "Basename for filtering input files. Only files starting with this are used. It is also prepended to result files (optional)")
            ("sourcepath,s", program_options::value<string>()->default_value("."), "Current cluster file (netCDF)")
            ("exclude-degenerates,n", program_options::value<bool>()->default_value(true), "Exclude results of tracks of length one")
            ("create-length-statistics,1", "Create a statistic of track lengths.")
            ("length-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(length_hist_default), "List of track-length values for histogram bins")
            ("create-speed-statistics,2", "Evaluate speeds of clusters in tracks, based on geometric center point")
            ("speed-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(speed_hist_default), "Speed histogram. Values in [m/s]")
            ("create-direction-statistics,3", "Evaluate directions of clusters in tracks, based on geometric center point")
            ("direction-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(direction_hist_default), "Direction histogram. Values in [deg]. Use with radolan grid only!!")
            ("create-cluster-statistics,4", "Evaluate each cluster in each track in terms of size.")
            ("cluster-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(cluster_hist_default), "List of cluster size values for histogram bins")
            ("create-cumulated-size-statistics,5", "Evaluate each track in terms of cumulative size. Warning: this process takes a lot of memory.")
            ("size-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(size_hist_default), "List of cumulated track size values for histogram bins")
            ("write-track-dictionary,t", "Write out a dictionary listing tracks with number of clusters etc.")
#if WITH_VTK
            ("write-center-tracks-as-vtk,e", "Write tracks out as .vtk files")
            ("write-cumulated-tracks-as-vtk,m", "Write cumulated tracks out as .vtk files. Only has effect if --create-cumulated-size-statistics is used")
            ("vtk-dimensions", program_options::value<svec_t>()->multitoken(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
#endif    
            ("write-gnuplot-files,g", "write individual files for the statistics fit for use with gnuplot")
            ;

    program_options::variables_map vm;

    try
    {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    } catch (std::exception &e)
    {
        cerr << "FATAL:Error parsing command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help") == 1 || argc < 2)
    {
        cout << desc << "\n";
        return EXIT_SUCCESS;
    }

    if (vm.count("version") != 0)
    {
        cout << m3D::VERSION << endl;
        return EXIT_SUCCESS;
    }

    // Evaluate user input

    string basename, sourcepath;
    bool exclude_degenerates = true;

    bool create_length_statistics = false;
    bin_t length_histogram_classes;

    bool create_cluster_statistics = false;
    bin_t cluster_histogram_classes;

    bool create_speed_statistics = false;
    fvec_t speed_histogram_classes;

    bool create_direction_statistics = false;
    fvec_t direction_histogram_classes;

    bool create_cumulated_size_statistics = false;
    bin_t size_histogram_classes;

    bool write_center_tracks_as_vtk = false;
    bool write_cumulated_tracks_as_vtk = false;
    svec_t vtk_dim_names;

    bool write_gnuplot_files = false;

    bool write_track_dictionary = false;

    try
    {
        parse_commmandline(vm,
                basename,
                sourcepath,
                vtk_dim_names,
                create_length_statistics,
                length_histogram_classes,
                create_speed_statistics,
                speed_histogram_classes,
                create_direction_statistics,
                direction_histogram_classes,
                create_cluster_statistics,
                cluster_histogram_classes,
                exclude_degenerates,
                create_cumulated_size_statistics,
                size_histogram_classes,
                write_center_tracks_as_vtk,
                write_cumulated_tracks_as_vtk,
                write_gnuplot_files,
                write_track_dictionary);

    } catch (const std::exception &e)
    {
        cerr << "FATAL:" << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // initialize values

    size_t number_of_degenerates = 0;
    bin_t length_histogram(length_histogram_classes.size(), 0);
    bin_t size_histogram(size_histogram_classes.size(), 0);
    bin_t cluster_histogram(cluster_histogram_classes.size(), 0);
    bin_t speed_histogram(speed_histogram_classes.size(), 0);
    bin_t direction_histogram(direction_histogram_classes.size(), 0);

    // This map contains a mapping from the id type to lists of clusters.

    Track<FS_TYPE>::trackmap track_map;

    // Collect all files at sourcepath that start with basename and
    // end with -clusters.nc

    namespace fs = boost::filesystem;

    fs::path source_path(sourcepath);

    size_t spatial_rank = 0;
    size_t value_rank = 0;

    std::string coords_filename;

    cout << "Collecting track data from NetCDF files: " << endl;
    
    // This file should not be closed until the code has
    // run through, or netCDF will upchuck exceptions when
    // accessing coordinate system
    netCDF::NcFile *coords_file = NULL;
    vector<NcDim> dimensions;
    vector<NcVar> dimension_vars;
    vector<string> dim_names;
    vector<string> var_names;

    size_t c_id = 0;

    if (fs::is_directory(source_path))
    {
        fs::directory_iterator dir_iter(source_path);
        fs::directory_iterator end;

        // Iterate over the "-clusters.nc" - files in sourcepath
        
        while (dir_iter != end)
        {
            double time_reading = 0;
            double time_searching = 0;
            double time_inserting = 0;
            double time_pushing_clusters = 0;
            double time_deleting = 0;
            double time_constructing = 0;

            fs::path f = dir_iter->path();

            if (fs::is_regular_file(f) && boost::algorithm::ends_with(f.filename().generic_string(), "-clusters.nc"))
            {
                // read the ClusterList from the file

#if APPLE
                double free_mem, used_mem;
                get_memory(free_mem, used_mem);
                cout << "Processing " << f.filename().generic_string()
                        << " ( used " << used_mem << " mb, "
                        << free_mem << " mb available) ... ";
#else 
#endif

                // Remember the first one

                if (coords_filename.empty())
                {
                    coords_filename = f.generic_string();
                    coords_file = new NcFile(coords_filename, NcFile::read);
                    
                    // Read "featurespace_dimensions"
                    cout << "Constructing coordinate system from file " << coords_filename << endl;

                    string fs_dimensions;
                    coords_file->getAtt("featurespace_dimensions").getValues(fs_dimensions);
                    dim_names = vectors::from_string<string>(fs_dimensions);
                    
                    string fs_vars;
                    coords_file->getAtt("featurespace_variables").getValues(fs_vars);

                    vector<string> all_var_names = vectors::from_string<string>(fs_vars);
                    for (size_t i = spatial_rank; i < all_var_names.size(); i++)
                    {
                        var_names.push_back(all_var_names[i]);
                    }

                    #if WITH_VTK
                    VisitUtils<FS_TYPE>::update_vtk_dimension_mapping(dim_names, vtk_dim_names);
                    #endif
                }

                try
                {
                    // Read the cluster list from file

                    // TODO: check if the actual dimensions and their 
                    // order remain constant

                    start_timer();
                    ClusterList<FS_TYPE>::ptr cluster_list = ClusterList<FS_TYPE>::read(f.generic_string());
                    time_reading += stop_timer();
                    
                    cout << "Processing " << f.filename().generic_string()
                        << " (" << cluster_list->size() << " clusters) ... ";

                    if (spatial_rank == 0) {
                        spatial_rank = cluster_list->dimensions.size();
                    } else if (spatial_rank != cluster_list->dimensions.size()) {
                        cerr << "FATAL:spatial range must remain identical across the track" << endl;
                        return EXIT_FAILURE;
                    }

                    // TODO: check if the actual variables and their order
                    // remain constant

                    size_t v_rank = cluster_list->feature_variables.size() - spatial_rank;
                    if (value_rank == 0) {
                        value_rank = v_rank;
                    } else if (v_rank != value_rank) {
                        cerr << "FATAL:value range must remain identical across the track" << endl;
                        return EXIT_FAILURE;
                    }

                    // Iterate over the clusters in the list we just read
                    Cluster<FS_TYPE>::list::const_iterator ci;
                    for (ci = cluster_list->clusters.begin(); ci != cluster_list->clusters.end(); ++ci)
                    {
                        Cluster<FS_TYPE>::ptr cluster = (*ci);
                        m3D::id_t id = cluster->id;
                        Track<FS_TYPE>::ptr tm = NULL;
                        Track<FS_TYPE>::trackmap::const_iterator ti;

                        start_timer();
                        ti = track_map.find(id);
                        time_searching += stop_timer();
                            
                        if (ti == track_map.end()) {
                            // new entry
                            tm = new Track<FS_TYPE>();
                            tm->id = id;
                            start_timer();
                            track_map[id] = tm;
                            time_inserting += stop_timer();
                        } else {
                            tm = ti->second;
                        }

                        boost::filesystem::path sf(cluster_list->source_file);
                        tm->sourcefiles.push_back(sf.filename().generic_string());

                        // Instead of the original cluster, use a TrackCluster, which
                        // has facilities of writing it's point list to disk and read
                        // it back on demand, saving memory.

                        start_timer();
                        TrackCluster<FS_TYPE>::ptr tc = new TrackCluster<FS_TYPE>(c_id++, cluster);
                        tm->clusters.push_back(tc);
                        time_pushing_clusters += stop_timer();

                        // dispose of the points in the original cluster
                        start_timer();
                        cluster->clear(true);
                        time_deleting += stop_timer();
                    }
                    
                    start_timer();
                    delete cluster_list;
                    time_deleting += stop_timer();
                    
                    cout << "tracking map has " << track_map.size() << " tracks" << endl;

                } catch (netCDF::exceptions::NcException &e)
                {
                    cerr << "ERROR:" << e.what() << endl;
                } catch (std::exception &e)
                {
                    cerr << "ERROR:" << e.what() << endl;
                }

                cout << "done." << endl;
            }
            
            cout << "Time spend:" << endl;
            cout << "\t reading:" << time_reading << endl;
            cout << "\t searching:" << time_searching << endl;
            cout << "\t inserting:" << time_inserting << endl;
            cout << "\t pushing clusters:" << time_pushing_clusters << endl;
            cout << "\t deleting:" << time_deleting << endl;

            dir_iter++;
        }

        cout << endl;

        cout << endl;
        cout << "Collating dimensions and dimension data ..." << endl;

        // Construct coordinate system

        CoordinateSystem<FS_TYPE> *coord_system = NULL;

        cout << "Spatial range variables (dimensions): " << dim_names << endl;
        cout << "Value range variables: " << var_names << endl;

        multimap<string, NcVar> vars = coords_file->getVars();
        multimap<string, NcVar>::iterator vi;

        // Only use dimensions, that have a variable of the
        // exact same name

        for (size_t di = 0; di < dim_names.size(); di++)
        {
            NcDim dim = coords_file->getDim(dim_names[di]);

            if (dim.isNull())
            {
                cerr << "FATAL: dimension " << dim_names[di] << " does not exist " << endl;
                exit(EXIT_FAILURE);
            }

            NcVar var = coords_file->getVar(dim_names[di]);

            if (var.isNull())
            {
                cerr << "FATAL: variable " << dim_names[di] << " does not exist " << endl;
                exit(EXIT_FAILURE);
            }

            dimensions.push_back(dim);
            dimension_vars.push_back(var);
        }

        coord_system = new CoordinateSystem<FS_TYPE>(dimensions, dimension_vars);

        // Now we have a cluster map. Let's print it for debug purposes

        cout << endl;
        cout << "Keying up tracks: " << endl;

        Track<FS_TYPE>::trackmap::iterator tmi;

        for (tmi = track_map.begin(); tmi != track_map.end(); ++tmi)
        {
            Track<FS_TYPE>::ptr track = tmi->second;

            cout << "Track #" << tmi->first
                    << " (" << track->clusters.size() << " clusters)"
                    << endl;

            size_t i = 0;

            std::list<Cluster<FS_TYPE>::ptr>::const_iterator ti;

            for (ti = track->clusters.begin(); ti != track->clusters.end(); ++ti)
            {
                Cluster<FS_TYPE>::ptr c = *ti;
                cout << "  [" << i++ << "] x=" << c->geometrical_center(spatial_rank) << endl;
                c->clear(true);
            }
        }

#if WITH_VTK
        if (write_center_tracks_as_vtk)
        {
            // Write center tracks
            VisitUtils<FS_TYPE>::write_center_tracks_vtk(track_map,
                    basename,
                    spatial_rank,
                    exclude_degenerates);
        }
#endif        
        // dictionary?

        if (write_track_dictionary)
        {
            ofstream dict("track-dictionary.json");

            dict << "{" << endl;
            dict << "  \"spatial_range\":" << to_json(dim_names) << "," << endl;
            dict << "  \"value_range\":" << to_json(var_names) << "," << endl;
            dict << "  \"number_of_tracks\":" << track_map.size() << "," << endl;

            dict << "  " << "\"tracks\":[" << endl;

            Track<FS_TYPE>::trackmap::iterator tmi;

            size_t track_index = 0;

            for (tmi = track_map.begin(); tmi != track_map.end(); ++tmi)
            {
                Track<FS_TYPE>::ptr track = tmi->second;

                dict << "    {" << endl;

                // Identifier
                dict << "      \"id\":" << track->id << "," << endl;

                // Length
                dict << "      \"length\":" << track->clusters.size() << "," << endl;

                // Clusters
                dict << "      \"clusters\":[" << endl;

                track->min.resize(value_rank, std::numeric_limits<FS_TYPE>::max());
                track->max.resize(value_rank, std::numeric_limits<FS_TYPE>::min());

                size_t i=0;
                std::list<Cluster<FS_TYPE>::ptr>::const_iterator ti;
                for (ti = track->clusters.begin(); ti != track->clusters.end(); ++ti) {
                    Cluster<FS_TYPE>::ptr c = (*ti);

                    vector<FS_TYPE> min, max, median;
                    c->variable_ranges(min, max, median);

                    dict << "        {" << endl;
                    dict << "          \"size\":" << c->size() << "," << endl;
                    dict << "          \"sourcefile\":\"" << track->sourcefiles[i] << "\"," << endl;
                    dict << "          \"geometrical_center\":" << to_json(c->geometrical_center(spatial_rank)) << "," << endl;
                    dict << "          \"min\":" << to_json(min) << "," << endl;
                    dict << "          \"max\":" << to_json(max) << "," << endl;
                    dict << "          \"median\":" << to_json(median) << endl;
                    dict << "          \"has_margin_points\":" << (c->has_margin_points() ? "Y" : "N") << endl;
                    dict << "        }";

                    // min/max
                    for (size_t k = 0; k < value_rank; k++)
                    {
                        if (min[k] < track->min[k]) track->min[k] = min[k];
                        if (max[k] > track->max[k]) track->max[k] = max[k];
                    }

                    if (i < (track->clusters.size() - 1))
                        dict << "," << endl;

                    // Dispose of the points again to free up memory.
                    c->clear(true);
                    
                    i++;
                }

                dict << "      ]," << endl;

                // Limits / Median
                dict << "      \"min\":" << to_json(track->min) << "," << endl;
                dict << "      \"max\":" << to_json(track->max) << endl;
                dict << "    }";

                if (track_index < (track_map.size() - 1))
                {
                    dict << "," << endl;
                }

                track_index++;
            }

            dict << "  ]";
            dict << "}" << endl;
            dict.close();
        }


        // Write cumulated tracks as netcdf and vtk files

        cout << "Cumulating track data: " << endl;

        map<size_t, size_t> track_lengths;
        map<size_t, size_t> track_sizes;
        map<size_t, size_t> cluster_sizes;
        map<float, size_t> speeds;
        map<float, size_t> directions;

        // Iterate over the collated tracks

        for (tmi = track_map.begin(); tmi != track_map.end(); ++tmi)
        {
            size_t points_processed = 0;

            Track<FS_TYPE>::ptr track = tmi->second;

            size_t track_length = track->clusters.size();

            if (track_length == 1)
            {
                number_of_degenerates++;
            }

            if (exclude_degenerates && track_length == 1)
            {
                continue;
            }

            cout << "Processing track #" << tmi->first << " (" << track->clusters.size() << " clusters)" << endl;

            // count up histogram

            if (create_length_statistics)
            {
                bool exceeded_max_class = false;

                add_value_to_histogram(length_histogram_classes, length_histogram, track_length, exceeded_max_class);

                // add to distribution

                map<size_t, size_t>::iterator tlfi = track_lengths.find(track_length);
                if (tlfi == track_lengths.end())
                    track_lengths[track_length] = 1;
                else
                    track_lengths[track_length] = (tlfi->second + 1);
            }

            // Skip the rest if cumulative size statistics and cluster / speed
            // stats are all off

            if (!create_cumulated_size_statistics
                    && !create_cluster_statistics
                    && !create_speed_statistics
                    && !create_direction_statistics)
                continue;

            // For each track, create an array index. Start with empty index.

            ArrayIndex<FS_TYPE> index(coord_system->get_dimension_sizes(), false);

            // Iterate over the clusters in the track and sum up

            std::list<Cluster<FS_TYPE>::ptr>::iterator ti;
            Cluster<FS_TYPE>::ptr previous_cluster = NULL;

            for (ti = track->clusters.begin(); ti != track->clusters.end(); ++ti)
            {
                Cluster<FS_TYPE>::ptr cluster = *ti;

                if (create_cluster_statistics)
                {
                    size_t cluster_size = cluster->size();

                    // add to histogram

                    bool exceeded_max_class = false;

                    add_value_to_histogram(cluster_histogram_classes, cluster_histogram, cluster_size, exceeded_max_class);

                    if (exceeded_max_class)
                    {
                        cout << "Cluster #" << cluster->id
                                << " (size " << cluster->size() << ")"
                                << " exceeded max cluster size"
                                << cluster_histogram_classes.back() << endl;
                    }

                    // add to distribution

                    map<size_t, size_t>::iterator csfi = cluster_sizes.find(cluster_size);

                    if (csfi == cluster_sizes.end())
                    {
                        cluster_sizes[cluster_size] = 1;
                    } else
                    {
                        cluster_sizes[cluster_size] = (csfi->second + 1);
                    }
                }

                if (create_speed_statistics)
                {
                    if (previous_cluster != NULL)
                    {
                        // calculate speed in m/s. All vectors 2D at this time

                        vector<FS_TYPE> p1 = previous_cluster->geometrical_center(2);
                        vector<FS_TYPE> p2 = cluster->geometrical_center(2);

                        // RADOLAN is in km. -> Tranform to meters
                        FS_TYPE dS = vector_norm(p2 - p1) * 1000;

                        // TODO: this information must come from somewhere!!
                        // => Cluster files need dT info! See ticket #227

                        FS_TYPE dT = 300.0;

                        // Average speed
                        FS_TYPE speed = dS / dT;

                        bool exceeded_max_class = false;

                        add_value_to_histogram<float>(speed_histogram_classes, speed_histogram, speed, exceeded_max_class);

                        if (exceeded_max_class)
                        {
                            cout << "Cluster #" << cluster->id
                                    << " (speed " << speed << " m/s)"
                                    << " exceeded max speed "
                                    << speed_histogram_classes.back() << " m/s"
                                    << endl;
                        }

                        map<float, size_t>::iterator si = speeds.find(track_length);

                        if (si == speeds.end())
                        {
                            speeds[speed] = 1;
                        } else
                        {
                            speeds[speed] = (si->second + 1);
                        }
                    }
                }

                // Directional Statistics
                if (create_direction_statistics)
                {
                    // Using the standard nationwide composite here
                    // NOTE: by doing this, the code is made dependant of the
                    // radolan code and fixed to the radolan coordinate system
                    // (which is BAD).

                    // TODO: find a more generic way to deal with this problem

                    RDCoordinateSystem *rcs = new RDCoordinateSystem(RD_RX);

                    if (previous_cluster != NULL)
                    {
                        // calculate speed in m/s. All vectors 2D at this time

                        vector<FS_TYPE> p1 = previous_cluster->geometrical_center(2);
                        vector<FS_TYPE> p2 = cluster->geometrical_center(2);

                        vector<FS_TYPE> dP = p2 - p1;

                        RDCartesianPoint p;

                        // Assuming --vtk-dimensions=x,y

#if WITH_VTK
                        if (!::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES.empty())
                        {
                            p.x = p1.at(::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES.at(0));
                            p.y = p1.at(::m3D::utils::VisitUtils<FS_TYPE>::VTK_DIMENSION_INDEXES.at(1));
                        } else
                        {
#endif
                            // traditionally the coordinates in NetCDF files
                            // come in the z,y,x order. Thus pick the default
                            // indexes accordingly

                            p.x = p1.at(1);
                            p.y = p1.at(0);
#if WITH_VTK
                        }
#endif                        
                        // Obtain geographical coordinate
                        RDGeographicalPoint c = rcs->geographicalCoordinate(p);

                        // move north 1 degree
                        RDGeographicalPoint cn = c;
                        cn.latitude = cn.latitude + 1.0;

                        // Transform back
                        RDCartesianPoint pn = rcs->cartesianCoordinate(cn);

                        // obtain difference vector and normalize it.
                        // this will point north

                        vector<FS_TYPE> n(2);
                        n[0] = pn.x - p.x;
                        n[1] = pn.y - p.y;

                        n = n / vector_norm(n);

                        // obtain the angle beta between the radolan grid
                        // and north at the distance vector's origin

                        vector<FS_TYPE> ey(2);
                        ey[0] = 0;
                        ey[1] = 1;

                        float beta = acos(ey * n);

                        // calculate angle alpha between the distance vector
                        // and the grid

                        float alpha = acos((ey * dP) / vector_norm(dP));

                        // the direction with north is the sum of the two angles
                        // in DEG

                        float direction = 180.0 * (alpha + beta) / M_PI_2;

                        if (direction > 360.0)
                        {
                            direction -= 360.0;
                        }

                        bool exceeded_max_class = false;

                        add_value_to_histogram<float>(direction_histogram_classes, direction_histogram, direction, exceeded_max_class);

                        map<float, size_t>::iterator si = directions.find(direction);

                        if (si == directions.end())
                        {
                            directions[direction] = 1;
                        } else
                        {
                            directions[direction] = (si->second + 1);
                        }
                    }
                }

                previous_cluster = cluster;

                // Skip the rest if cumulative size statistics are off

                if (!create_cumulated_size_statistics)
                    continue;

                // Iterate over the points of the cluster and add points
                // to the array index. If the point exists, add it's values
                // to the existing one. If not, add a new point. In the end
                // the index contains all points that all clusters in the
                // track occupied at any time

                Point<FS_TYPE>::list::iterator pi;

                for (pi = cluster->get_points().begin(); pi != cluster->get_points().end(); ++pi)
                {
                    Point<FS_TYPE>::ptr p = *pi;

                    Point<FS_TYPE>::ptr indexed = index.get(p->gridpoint);

                    if (indexed == NULL)
                    {
                        // this makes a copy of the point in the index
                        index.set(p->gridpoint, p);

                        // get the copied point
                        indexed = index.get(p->gridpoint);
                    }

                    // only add up the value range

                    for (size_t k = 0; k < value_rank; k++)
                    {
                        indexed->values[spatial_rank + k] += p->values[spatial_rank + k];
                    }

                    points_processed++;
                }

                // Purge points in memory
                cluster->clear(true);
            }

            // A track's size is the number of non-null points in
            // the index

            size_t track_size = index.count();

#if WITH_VTK
            if (write_cumulated_tracks_as_vtk)
            {
                Point<FS_TYPE>::list cumulatedList;
                index.replace_points(cumulatedList);
                string vtk_path = basename + "_cumulated_track_" + boost::lexical_cast<string>(tmi->first) + ".vtk";
                VisitUtils<FS_TYPE>::write_pointlist_all_vars_vtk(vtk_path, &cumulatedList, vector<string>());
            }
#endif
            // Delete the index to save memory
            index.clear(true);

            // time to delete the map

            map<size_t, size_t>::iterator tlfi = track_sizes.find(track_size);

            if (tlfi == track_sizes.end())
            {
                track_sizes[track_size] = 1;
            } else
            {
                track_sizes[track_size] = (tlfi->second + 1);
            }

            // Update histogram

            bool exceeded_max_class = false;

            add_value_to_histogram(size_histogram_classes, size_histogram, track_size, exceeded_max_class);

            // Write out the cumulative file as netcdf and vtk

            // Don't need it anymore

            delete track;

            // NETCDF

            cout << "  (processed " << points_processed << " points)" << endl;
        }

        cout << "done." << endl;

        // Write out statistical data in file file(s)

        string file_fn = basename + "_trackstats.txt";

        ofstream file(file_fn.c_str());

        // number of tracks

        file << "Tracking report for basename " + basename << endl;
        cout << "Tracking report for basename " + basename << endl;

        file << "(degenerates are excluded: " << (exclude_degenerates ? "yes" : "no") << ")" << endl;
        cout << "(degenerates are excluded: " << (exclude_degenerates ? "yes" : "no") << ")" << endl;

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

        if (create_length_statistics)
        {
            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            file << "Track length distribution" << endl;
            cout << "Track length distribution" << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            double avg = average<size_t>(track_lengths, true);
            file << "(Average track length = " << avg << ")" << endl;
            cout << "(Average track length = " << avg << ")" << endl;

            map<size_t, size_t>::reverse_iterator rend = track_lengths.rbegin();
            file << "(Maximum track length = " << rend->first << ")" << endl;
            cout << "(Maximum track length = " << rend->first << ")" << endl;

            file << "length,number" << endl;
            cout << "length,number" << endl;

            map<size_t, size_t>::iterator si;

            for (si = track_lengths.begin(); si != track_lengths.end(); si++)
            {
                file << si->first << "," << si->second << endl;
                cout << si->first << "," << si->second << endl;
            }

            file << endl;

            // Histogram

            file << endl;
            cout << endl;

            file << "Length Histogram:" << endl;
            cout << "Length Histogram:" << endl;

            file << "max length,number" << endl;
            cout << "max length,number" << endl;

            print_histogram(length_histogram_classes, length_histogram, file);

            if (write_gnuplot_files)
            {
                write_histogram("lengths-hist.txt", "number of tracks", "track length", length_histogram_classes, length_histogram);
                write_values<size_t>("lengths.txt", "number of tracks", "track length", track_lengths);
            }
        }

        //
        // speed
        //

        if (create_speed_statistics)
        {
            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            file << "Speed distribution" << endl;
            cout << "Speed distribution" << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            double avg = average<float>(speeds, true);
            file << "(Average speed = " << avg << ")" << endl;
            cout << "(Average speed = " << avg << ")" << endl;

            map<float, size_t>::reverse_iterator srend = speeds.rbegin();
            file << "(Maximum speed = " << srend->first << ")" << endl;
            cout << "(Maximum speed = " << srend->first << ")" << endl;

            file << "speed,number" << endl;
            cout << "speed,number" << endl;

            map<float, size_t>::iterator spi;

            for (spi = speeds.begin(); spi != speeds.end(); spi++)
            {
                file << spi->first << "," << spi->second << endl;
                cout << spi->first << "," << spi->second << endl;
            }

            file << endl;

            // Histogram

            file << endl;
            cout << endl;

            file << "Speed Histogram:" << endl;
            cout << "Speed Histogram:" << endl;

            file << "speed class,number" << endl;
            cout << "speed class,number" << endl;

            print_histogram<float>(speed_histogram_classes, speed_histogram, file);

            if (write_gnuplot_files)
            {
                write_histogram("speeds-hist.txt", "number of clusters", "speed [m/s]", speed_histogram_classes, speed_histogram);
                write_values<float>("speeds.txt", "number of clusters", "speed in [m/s]", speeds);
            }

        }

        //
        // directions
        //

        if (create_direction_statistics)
        {
            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            file << "Direction distribution" << endl;
            cout << "Direction distribution" << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            double avg = average<float>(directions, true);
            file << "(Average direction = " << avg << ")" << endl;
            cout << "(Average direction = " << avg << ")" << endl;

            file << "direction,number" << endl;
            cout << "direction,number" << endl;

            map<float, size_t>::iterator spi;

            for (spi = directions.begin(); spi != directions.end(); spi++)
            {
                file << spi->first << "," << spi->second << endl;
                cout << spi->first << "," << spi->second << endl;
            }

            file << endl;

            // Histogram

            file << endl;
            cout << endl;

            file << "Direction Histogram:" << endl;
            cout << "Direction Histogram:" << endl;

            file << "direction class,number" << endl;
            cout << "direction class,number" << endl;

            print_histogram<float>(direction_histogram_classes, direction_histogram, file);

            if (write_gnuplot_files)
            {
                write_histogram("directions-hist.txt", "number of clusters", "direction in [deg]", direction_histogram_classes, direction_histogram);
                write_values<float>("directions.txt", "number of clusters", "direction in [deg]", directions);
            }

        }

        //
        // Cluster sizes
        //

        if (create_cluster_statistics)
        {
            file << endl;
            cout << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            file << "Cluster size distribution" << endl;
            cout << "Cluster size distribution" << endl;

            file << "------------------------------------------------" << endl;
            cout << "------------------------------------------------" << endl;

            double avg = average<size_t>(cluster_sizes, true);
            file << "(Average cluster size = " << avg << ")" << endl;
            cout << "(Average cluster size = " << avg << ")" << endl;

            map<size_t, size_t>::reverse_iterator rend = cluster_sizes.rbegin();
            file << "(Maximum cluster size = " << rend->first << ")" << endl;
            cout << "(Maximum cluster size = " << rend->first << ")" << endl;

            file << "size,number" << endl;
            cout << "size,number" << endl;

            map<size_t, size_t>::iterator si;

            for (si = cluster_sizes.begin(); si != cluster_sizes.end(); si++)
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

            print_histogram(cluster_histogram_classes, cluster_histogram, file);

            if (write_gnuplot_files)
            {
                write_histogram("sizes-hist.txt", "number of clusters", "size [#gridpoints]", cluster_histogram_classes, cluster_histogram);
                write_values("sizes.txt", "number of clusters", "size [#gridpoints]", cluster_sizes);
            }

        }

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

            double avg = average<size_t>(track_sizes, false);
            file << "(Average track size = " << avg << ")" << endl;
            cout << "(Average track size = " << avg << ")" << endl;

            map<size_t, size_t>::reverse_iterator rend = track_sizes.rbegin();
            file << "(Maximum track size = " << rend->first << ")" << endl;
            cout << "(Maximum track size = " << rend->first << ")" << endl;

            file << "size,number" << endl;
            cout << "size,number" << endl;

            map<size_t, size_t>::iterator si;

            for (si = track_sizes.begin(); si != track_sizes.end(); si++)
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

            print_histogram(size_histogram_classes, size_histogram, file);

            if (write_gnuplot_files)
            {
                write_histogram("cumulative-sizes-hist.txt", "number of tracks", "size [#gridpoints]", size_histogram_classes, size_histogram);
                write_values<size_t>("cumulative-sizes.txt", "number of tracks", "size [#gridpoints]", track_sizes);
            }
        }

        // close the cluster file we used to get the coordinate
        // system data

        delete coords_file;
    } else
    {
        cerr << "Argument --sourcepath does not point to a directory (sourcepath=" + sourcepath + ")" << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
};
