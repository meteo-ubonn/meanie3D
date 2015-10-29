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

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <locale>
#include <limits>
#include <stdlib.h>
#include <netcdf>

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
        string &sourcepath,
        bool &create_length_stats,
        bin_t &length_histogram_bins,
        bool &create_speed_stats,
        fvec_t &speed_histogram_bins,
        bool & create_direction_stats,
        fvec_t &direction_histogram_classes,
        bool &create_cluster_stats,
        bin_t &cluster_histogram_bins,
        bool &exclude_degenerates,
        bool &write_center_tracks_as_vtk,
        bool &write_gnuplot_files,
        bool &write_track_dictionary)
{
    sourcepath = vm["filename"].as<string>();

    write_center_tracks_as_vtk = vm.count("write-center-tracks-as-vtk") > 0;

    write_gnuplot_files = vm.count("write-gnuplot-files") > 0;

    write_track_dictionary = vm.count("write-track-dictionary") > 0;

    exclude_degenerates = vm.count("exclude-degenerates") > 0;

    create_length_stats = vm.count("create-length-statistics") > 0;
    length_histogram_bins = vm["length-histogram-classes"].as<bin_t>();

    create_speed_stats = vm.count("create-speed-statistics") > 0;
    speed_histogram_bins = vm["speed-histogram-classes"].as<fvec_t>();

    create_direction_stats = vm.count("create-direction-statistics") > 0;
    direction_histogram_classes = vm["direction-histogram-classes"].as<fvec_t>();

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
void add_value_to_histogram(const vector<T> &classes, bin_t &counts, T value, bool &exceeded_max_class)
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
void print_histogram(const vector<T> &classes, const bin_t &values, ofstream &file)
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

    file << "#" << x << "\t" << y << endl;

    for (size_t i = 0; i < classes.size(); i++)
    {
        file << classes[i] << "\t" << values[i] << endl;
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

    file << "track length\tnumber of tracks" << endl;

    typename map<T, size_t>::const_iterator si;

    for (si = values.begin(); si != values.end(); si++)
    {
        file << si->first << "\t" << si->second << endl;
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
#pragma mark main

/**
 *
 *
 */
int main(int argc, char** argv)
{
    using namespace m3D;
    using namespace m3D::utils;
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
            ("filename,f", program_options::value<string>()->default_value("."), "Cluster file (KONRAD format)")
            ("exclude-degenerates", "Exclude results of tracks of length one")
            ("create-length-statistics", "Create a statistic of track lengths.")
            ("length-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(length_hist_default), "List of track-length values for histogram bins")
            ("create-speed-statistics", "Evaluate speeds of clusters in tracks, based on geometric center point")
            ("speed-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(speed_hist_default), "Speed histogram. Values in [m/s]")
            ("create-direction-statistics", "Evaluate directions of clusters in tracks, based on geometric center point")
            ("direction-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(direction_hist_default), "Direction histogram. Values in [deg]. Use with radolan grid only!!")
            ("create-cluster-statistics", "Evaluate each cluster in each track in terms of size.")
            ("cluster-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(cluster_hist_default), "List of cluster size values for histogram bins")
#if WITH_VTK
            ("write-center-tracks-as-vtk", "Write tracks out as .vtk files")
#endif
            ("write-track-dictionary", "Write out a dictionary listing tracks with number of clusters etc.")
            ("write-gnuplot-files,g", "write individual files for the statistics fit for use with gnuplot")
            ;

    program_options::variables_map vm;

    try
    {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    }    catch (std::exception &e)
    {
        cerr << "FATAL:Error parsing command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help") == 1 || argc < 2)
    {
        cout << desc << "\n";
        return 1;
    }

    // Version

    if (vm.count("version") != 0)
    {
        cout << m3D::VERSION << endl;
        exit(EXIT_FAILURE);
        ;
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

    //    bool create_cumulated_size_statistics = false;
    //    bin_t size_histogram_classes;

    bool write_center_tracks_as_vtk = false;
    //    bool write_cumulated_tracks_as_vtk = false;
    //    svec_t vtk_dim_names;

    bool write_gnuplot_files = false;

    bool write_track_dictionary = false;

    try
    {
        parse_commmandline(vm,
                sourcepath,
                create_length_statistics,
                length_histogram_classes,
                create_speed_statistics,
                speed_histogram_classes,
                create_direction_statistics,
                direction_histogram_classes,
                create_cluster_statistics,
                cluster_histogram_classes,
                exclude_degenerates,
                write_center_tracks_as_vtk,
                write_gnuplot_files,
                write_track_dictionary);

    }    catch (const std::exception &e)
    {
        cerr << "FATAL:" << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // initialize values

    size_t number_of_degenerates = 0;
    bin_t length_histogram(length_histogram_classes.size(), 0);
    //    bin_t size_histogram(size_histogram_classes.size(),0);
    bin_t cluster_histogram(cluster_histogram_classes.size(), 0);
    bin_t speed_histogram(speed_histogram_classes.size(), 0);
    bin_t direction_histogram(direction_histogram_classes.size(), 0);

    // This map contains a mapping from the id type to lists of clusters.

    typedef std::vector< ConradCluster<FS_TYPE> > track_t;
    typedef std::map<m3D::id_t, track_t* > trackmap_t;

    trackmap_t track_map;

    // Collect all files at sourcepath that start with basename and
    // end with -clusters.nc

    namespace fs = boost::filesystem;

    fs::path source_path(sourcepath);

    basename = fs::basename(sourcepath);

    cout << "Collecting track data from NetCDF files: " << endl;

    if (fs::exists(sourcepath))
    {
        // Re-using track_t here. Technically identical,
        // but don't get confused. This is not a track but
        // the list of all clusters from that file.

        track_t clusters = ConradCluster<FS_TYPE>::read_conrad_short(sourcepath);

        if (clusters.empty())
        {
            cout << "No clusters to process." << endl;
            exit(0);
        } else
        {
            cout << "Processing " << clusters.size() << " clusters" << endl;
        }

        track_t::const_iterator ci;

        for (ci = clusters.begin(); ci != clusters.end(); ++ci)
        {
            ConradCluster<FS_TYPE> cluster = (*ci);

            m3D::id_t id = cluster.id;

            trackmap_t::const_iterator ti = track_map.find(id);

            track_t *tm = NULL;

            if (ti == track_map.end())
            {
                // new entry

                tm = new track_t();

                track_map[id] = tm;
            } else
            {
                tm = ti->second;
            }

            // append the cluster to the end of the track. Note that this
            // is a copy operation. It's expensive, but it saves us from
            // having to keep track of open files and stuff.

            tm->push_back(cluster);
        }

        cout << endl;

        // Now we have a cluster map. Let's print it for debug purposes

        cout << "Keying up tracks: " << endl;

        for (trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            track_t *track = tmi->second;

            cout << "Track #" << tmi->first << " (" << track->size() << " clusters)" << endl;

            size_t i = 0;

            for (track_t::iterator ti = track->begin(); ti != track->end(); ti++)
            {
                cout << "\t[" << i++ << "] x=" << ti->center() << endl;
            }
        }

        //        // Write center tracks
#if WITH_VTK
        if (write_center_tracks_as_vtk)
        {
            ::m3D::utils::VisitUtils<FS_TYPE>::write_center_tracks_vtk(track_map, basename, exclude_degenerates);
        }
#endif
        // dictionary?

        if (write_track_dictionary)
        {
            ofstream dict("track-dictionary.txt");

            for (trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
            {
                track_t *track = tmi->second;

                dict << "Track #" << tmi->first << " (" << track->size() << " clusters)" << endl;

                size_t i = 0;

                for (track_t::iterator ti = track->begin(); ti != track->end(); ti++)
                {
                    dict << "\t[" << i++ << "] x=" << ti->center() << endl;
                }
            }

            dict.close();
        }


        // Write cumulated tracks as netcdf and vtk files

        cout << "Cumulating track data: " << endl;

        // track length

        map<size_t, size_t> track_lengths;

        map<size_t, size_t> track_sizes;

        map<size_t, size_t> cluster_sizes;

        map<float, size_t> speeds;

        map<float, size_t> directions;

        // Iterate over the collated tracks

        for (trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            track_t *track = tmi->second;

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

            if (create_length_statistics)
            {
                bool exceeded_max_class = false;

                add_value_to_histogram(length_histogram_classes, length_histogram, track_length, exceeded_max_class);

                // add to distribution

                map<size_t, size_t>::iterator tlfi = track_lengths.find(track_length);

                if (tlfi == track_lengths.end())
                {
                    track_lengths[track_length] = 1;
                } else
                {
                    track_lengths[track_length] = (tlfi->second + 1);
                }
            }

            // Skip the rest if cumulative size statistics and cluster / speed
            // stats are all off

            if (!create_cluster_statistics && !create_speed_statistics && !create_direction_statistics)
            {
                continue;
            }

            // For each track, create an array index. Start with empty index.

            // ArrayIndex<FS_TYPE> index(coord_system,false);

            // Iterate over the clusters in the track and sum up

            track_t::iterator ti;

            ConradCluster<FS_TYPE> *previous_cluster = NULL;

            for (ti = track->begin(); ti != track->end(); ++ti)
            {
                ConradCluster<FS_TYPE> *cluster = &(*ti);

                if (create_cluster_statistics)
                {
                    size_t cluster_size = cluster->numCorePixels;

                    // add to histogram

                    bool exceeded_max_class = false;

                    add_value_to_histogram(cluster_histogram_classes, cluster_histogram, cluster_size, exceeded_max_class);

                    if (exceeded_max_class)
                    {
                        cout << "Cluster #" << cluster->id
                                << " (size " << cluster->numCorePixels << ")"
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

                        vector<FS_TYPE> p1 = previous_cluster->center();
                        vector<FS_TYPE> p2 = cluster->center();

                        // RADOLAN is in km. -> Tranform to meters
                        FS_TYPE dS = vector_norm(p2 - p1) * 1000;

                        // TODO: this information must come from somewhere!!
                        // => Cluster files need dT info! See ticket #227

                        FS_TYPE dT = (cluster->secondsSinceEpoch() - previous_cluster->secondsSinceEpoch());

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
                    // (see #

                    RDCoordinateSystem *rcs = new RDCoordinateSystem(RD_RX);

                    if (previous_cluster != NULL)
                    {
                        // calculate speed in m/s. All vectors 2D at this time

                        vector<FS_TYPE> p1 = previous_cluster->center();
                        vector<FS_TYPE> p2 = cluster->center();

                        vector<FS_TYPE> dP = p1 - p2;

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
                            p.x = 0;
                            p.y = 1;
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

                        // the direction with north is the sum of the two vectors
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

            }

            // Don't need it anymore

            delete track;
        }

        cout << "done." << endl;

        // Write out statistical data in file file(s)

        string file_fn = basename + "conrad_trackstats.txt";

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
                write_histogram("conrad_lengths-hist.txt", "number of tracks", "track length", length_histogram_classes, length_histogram);
                write_values<size_t>("conrad_lengths.txt", "number of tracks", "track length", track_lengths);
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
                write_histogram("conrad_speeds-hist.txt", "number of clusters", "speed [m/s]", speed_histogram_classes, speed_histogram);
                write_values<float>("conrad_speeds.txt", "number of clusters", "speed in [m/s]", speeds);
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
                write_histogram("conrad_directions-hist.txt", "number of clusters", "direction in [deg]", direction_histogram_classes, direction_histogram);
                write_values<float>("conrad_directions.txt", "number of clusters", "direction in [deg]", directions);
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
                write_histogram("conrad_sizes-hist.txt", "number of clusters", "size [#gridpoints]", cluster_histogram_classes, cluster_histogram);
                write_values("conrad_sizes.txt", "number of clusters", "size [#gridpoints]", cluster_sizes);
            }
        }
    } else
    {
        cerr << "Argument --sourcepath does not point to a directory (sourcepath=" + sourcepath + ")" << endl;
        exit(EXIT_FAILURE);
    }


    return 0;
};
