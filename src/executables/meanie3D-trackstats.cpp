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

#include <algorithm>
#include <exception>
#include <fstream>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <netcdf>
#include <unistd.h>
#include <set>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

#include <meanie3D/meanie3D.h>
#include <radolan/radolan.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;
using namespace Radolan;

namespace fs = boost::filesystem;

#pragma mark -
#pragma mark type definitions

/** Feature-space data type */
typedef double FS_TYPE;

typedef vector<float> fvec_t;
typedef vector<size_t> bin_t;
typedef vector<string> svec_t;

typedef enum
{
    Continue, Split, Merge
} link_type_t;

/**
 * Node in a track graph.
 */
typedef struct
{
    m3D::uuid_t uuid;
    m3D::id_t id;
    unsigned int step;
    unsigned int size;
} node_t;

/**
 * Link in a track graph.
 */
typedef struct
{
    m3D::uuid_t source; // source uuid
    m3D::uuid_t target; // target uuid
    m3D::id_t id; // target id 
    link_type_t type; // split, merge, continue?
} link_t;

// Store the cluster min/max and median values so that the
// points can be disposed of early on.
typedef map< m3D::id_t, vector<FS_TYPE> > val_map_t;

// Parameters

typedef struct
{
    std::string sourcepath;
    std::string basename;
    svec_t vtk_dim_names;
    bool create_length_stats;
    bin_t length_histogram_bins;
    bool create_speed_stats;
    fvec_t speed_histogram_bins;
    bool create_direction_stats;
    fvec_t direction_histogram_bins;
    bool create_cluster_stats;
    bin_t cluster_histogram_bins;
    bool exclude_degenerates;
    bool create_cumulated_size_stats;
    bin_t size_histogram_bins;
    bool create_cumulated_tracking_stats;
    bool write_center_tracks_as_vtk;
    bool write_cumulated_tracks_as_vtk;
    bool write_gnuplot_files;
    bool write_track_dictionary;
} parameter_t;

/**
 * Context object passed between various methods.
 */
typedef struct
{
    // Command line parameters
    parameter_t params;

    // Control flags and other properties
    bool need_points;

    // NetCDF properties
    size_t spatial_rank;
    size_t value_rank;
    netCDF::NcFile *coords_file;
    vector<NcDim> dimensions;
    vector<NcVar> dimension_vars;
    vector<string> dim_names;
    vector<string> var_names;

    // Context properties
    CoordinateSystem<FS_TYPE> *coord_system;

    Track<FS_TYPE>::trackmap track_map;

    vector<node_t> nodes;
    vector<link_t> links;
    unsigned int step;

    ClusterList<FS_TYPE>::ptr cluster_list;     // current processed cluster list
    Track<FS_TYPE>::ptr track;                  // current processed track
    ArrayIndex<FS_TYPE> *track_index;           // current processed track's index
    TrackCluster<FS_TYPE> *cluster;             // current processed cluster
    TrackCluster<FS_TYPE> *previous_cluster;    // last processed cluster
    size_t points_processed;

    bin_t length_histogram;
    bin_t size_histogram;
    bin_t cluster_histogram;
    bin_t speed_histogram;
    bin_t direction_histogram;

    size_t average_cluster_size;
    size_t number_of_clusters;
    size_t number_of_degenerates;

    val_map_t cluster_min;
    val_map_t cluster_max;
    val_map_t cluster_median;

    map<size_t, size_t> track_lengths;
    map<size_t, size_t> track_sizes;
    map<size_t, size_t> cluster_sizes;
    map<float, size_t> speeds;
    map<float, size_t> directions;

} trackstats_context_t;

/**
 * Prepare the stats context.
 */
trackstats_context_t initialiseContext(parameter_t params) {
    trackstats_context_t ctx;
    ctx.spatial_rank = 0;
    ctx.value_rank = 0;
    ctx.cluster_list = NULL;
    ctx.track = NULL;
    ctx.track_index = NULL;
    ctx.cluster = NULL;
    ctx.previous_cluster = NULL;
    ctx.points_processed = 0;
    ctx.coord_system = NULL;
    ctx.number_of_degenerates = 0;
    ctx.length_histogram = bin_t(params.length_histogram_bins.size(), 0);
    ctx.size_histogram = bin_t(ctx.params.size_histogram_bins.size(), 0);
    ctx.cluster_histogram = bin_t(params.cluster_histogram_bins.size(), 0);
    ctx.speed_histogram = bin_t(params.speed_histogram_bins.size(), 0);
    ctx.direction_histogram = bin_t(params.direction_histogram_bins.size(), 0);
    ctx.coords_file = NULL;
    ctx.need_points = params.create_cumulated_size_stats
                      || params.write_cumulated_tracks_as_vtk;
    ctx.step = 0;
    ctx.average_cluster_size = 0;
    ctx.number_of_clusters = 0;
    ctx.params = params;
    return ctx;
}

#pragma mark -
#pragma mark Comand line parsing

void parse_commmandline(program_options::variables_map vm, parameter_t &p) {
    p.sourcepath = vm["sourcepath"].as<string>();
    if (vm.count("basename") == 0) {
        fs::path sourcePath(p.sourcepath);
        p.basename = ::m3D::utils::common_component(sourcePath.generic_string(), ".nc");
    } else {
        p.basename = vm["basename"].as<string>();
    }
    p.write_center_tracks_as_vtk = vm.count("write-center-tracks-as-vtk") > 0;
    p.write_gnuplot_files = vm.count("write-gnuplot-files") > 0;
    p.write_cumulated_tracks_as_vtk = vm.count("write-cumulated-tracks-as-vtk") > 0;
    p.write_track_dictionary = vm.count("write-track-dictionary") > 0;
    p.exclude_degenerates = vm["exclude-degenerates"].as<bool>();
    // Default is the dimensions
    if (vm.count("vtk-dimensions") > 0) {
        // parse dimension list
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep(",");
        string str_value = vm["vtk-dimensions"].as<string>();
        tokenizer dim_tokens(str_value, sep);
        for (tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter) {
            string name = *tok_iter;
            p.vtk_dim_names.push_back(name);
        }
    }
    p.create_length_stats = vm.count("create-length-statistics") > 0;
    p.length_histogram_bins = vm["length-histogram-classes"].as<bin_t>();
    p.create_speed_stats = vm.count("create-speed-statistics") > 0;
    p.speed_histogram_bins = vm["speed-histogram-classes"].as<fvec_t>();
    p.create_direction_stats = vm.count("create-direction-statistics") > 0;
    p.direction_histogram_bins = vm["direction-histogram-classes"].as<fvec_t>();
    p.create_cumulated_size_stats = vm.count("create-cumulated-size-statistics") > 0;
    p.size_histogram_bins = vm["size-histogram-classes"].as<bin_t>();
    p.create_cluster_stats = vm.count("create-cluster-statistics") > 0;
    p.cluster_histogram_bins = vm["cluster-histogram-classes"].as<bin_t>();
    p.create_cumulated_tracking_stats = vm.count("create-cumulated-tracking-stats") > 0;
}

#pragma mark -
#pragma mark Histogram helpers

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
        bool &exceeded_max_class) {
    bool found_bin = false;
    exceeded_max_class = false;
    for (size_t i = 0; i < classes.size() && !found_bin; i++) {
        if (i > 0) {
            if (value <= classes[i] && value > classes[i - 1]) {
                counts[i] = counts[i] + 1;
                found_bin = true;
            }
        } else {
            if (value <= classes[0]) {
                counts[i] = counts[i] + 1;
                found_bin = true;
            }
        }
    }

    // not found? add to the last one
    if (!found_bin) {
        // add the last class with value 0 if required
        if (counts.size() == classes.size()) {
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
        ofstream &file) {
    for (size_t i = 0; i < classes.size(); i++) {
        file << classes[i] << "," << values[i] << endl;
        cout << classes[i] << "," << values[i] << endl;
    }
    // 'Over' bucket?
    if (values.size() > classes.size()) {
        file << " > " << classes[classes.size() - 1] << "," << values.back() << endl;
        cout << " > " << classes[classes.size() - 1] << "," << values.back() << endl;
    }
}

template <typename T>
void write_histogram(const std::string &filename,
        const std::string &x,
        const std::string &y,
        const vector<T> &classes,
        const bin_t &values) {
    ofstream file(filename.c_str());
    file << "#" << x << "  " << y << endl;
    for (size_t i = 0; i < classes.size(); i++) {
        file << classes[i] << "  " << values[i] << endl;
    }
    file.close();
}

template <typename T>
void write_values(const std::string &filename,
        const std::string &x,
        const std::string &y,
        const map<T, size_t> &values) {
    ofstream file(filename.c_str());
    file << "track length  number of tracks" << endl;
    typename map<T, size_t>::const_iterator si;
    for (si = values.begin(); si != values.end(); si++) {
        file << si->first << "  " << si->second << endl;
    }
    file.close();
}

template <typename T>
double average(const vector<FS_TYPE> &v) {
    // calculate the mean value
    double sum = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
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
double average(const map<T, size_t> &m, bool ignore_one = true) {
    // calculate the mean value
    double sum = 0.0;
    size_t total_count = 0;
    typename std::map<T, size_t>::const_iterator mi;
    for (mi = m.begin(); mi != m.end(); mi++) {
        T key = mi->first;
        size_t val = mi->second;
        if (key == 1 && ignore_one) {
            continue;
        }
        total_count += val;
        sum += boost::numeric_cast<double>(mi->first * mi->second);
    }
    return sum / boost::numeric_cast<double>(total_count);
}

#pragma mark -
#pragma mark graph construction

/**
 * Compare two nodes by uuid.
 * @param a
 * @param b
 * @return 
 */
bool compareNodesByUuid(const node_t &a, const node_t &b) {
    return (a.uuid >= b.uuid) ? true : false;
}

bool compareNodesByStep(const node_t &a, const node_t &b) {
    return (a.step < b.step) ? true : false;
}

/**
 * Finds all nodes with the given id in the node list and returns
 * them in ascending order of step
 * 
 * @param nodes
 * @param id
 * @param only include nodes in this step
 * @return 
 */
bool
findNode(const vector<node_t> &nodes, const m3D::id_t &id, const unsigned int &step, node_t &node) {
    bool found = false;
    for (size_t i = 0; i < nodes.size() && !found; i++) {
        node_t n = nodes[i];
        if (n.id == id && n.step == step) {
            node.id = n.id;
            node.step = n.step;
            node.uuid = n.uuid;
            found = true;
        }
    }
    return found;
}

void printNodes(const vector<node_t> &nodes) {
    for (size_t i = 0; i < nodes.size(); i++) {
        node_t node = nodes[i];
        cout << "step=" << node.step
                << " id=" << node.id
                << " uuid=" << node.uuid
                << endl;
    }
}

/**
 * Adds the link if no link with the same source->target and type exists.
 * 
 * @param ctx
 * @param source
 * @param target
 */
void
addUniqueLink(trackstats_context_t &ctx, node_t source, node_t target, link_type_t type) {
    // Figure out if there is a link for this already. This can
    // happen as the result of a split or link with continued id
    bool haveLink = false;
    for (int li = 0; li < ctx.links.size() && !haveLink; li++) {
        link_t link = ctx.links[li];
        haveLink = (link.source == source.uuid 
                && link.target == target.uuid
                && link.type == type);
    }
    if (!haveLink) {
        link_t link;
        link.source = source.uuid;
        link.target = target.uuid;
        link.type = type;
        switch(type) {
            case Split:
            case Continue:
                link.id = source.id;
                break;
            case Merge:
                link.id = target.id;
                break;
        }
        ctx.links.push_back(link);
    }
}

void
addGraphNode(trackstats_context_t &ctx, const node_t &node) {

    // Construct a new node for the tree data
    ctx.nodes.push_back(node);

    // check for split event
    id_map_t::const_iterator mi;
    for (mi = ctx.cluster_list->splits.begin(); mi != ctx.cluster_list->splits.end(); ++mi) {
        m3D::id_t sourceId = mi->first;
        node_t source;
        if (findNode(ctx.nodes, sourceId, node.step-1, source)) {
            id_set_t::const_iterator si;
            for (si = mi->second.begin(); si != mi->second.end(); si++) {
                m3D::id_t targetId = *si;
                if (targetId == node.id) {
                    addUniqueLink(ctx, source, node, Split);
                }
            }
        }
    }

    // check for merge events
    for (mi = ctx.cluster_list->merges.begin(); mi != ctx.cluster_list->merges.end(); ++mi) {
        if (mi->first == node.id) {
            id_set_t::const_iterator si;
            // Add all merged clusters at once, since those
            // should all have been added by the time this
            // is called (processed in the order of steps)
            for (si = mi->second.begin(); si != mi->second.end(); si++) {
                node_t source;
                m3D::id_t sourceId = (*si);
                if (findNode(ctx.nodes, sourceId, node.step-1, source)) {
                    addUniqueLink(ctx, source, node, Merge);
                }
            }
        }
    }

    // check for continuation event
    id_set_t::const_iterator ti = find(ctx.cluster_list->tracked_ids.begin(),
            ctx.cluster_list->tracked_ids.end(), node.id);
    if (ti != ctx.cluster_list->tracked_ids.end()) {
        node_t source;
        if (findNode(ctx.nodes, node.id, node.step-1, source)) {
            addUniqueLink(ctx, source, node, Continue);
        }
    }
}

#pragma mark -
#pragma mark Reading cluster data

/**
 * Constructs coordinate system, gets dimension- and variable names
 * and sets up vtk dimensions. Only done once.
 *
 * @param ctx
 * @param filename
 */
void getFeaturespaceInfo(trackstats_context_t &ctx, const std::string &filename) {
    if (ctx.coord_system == NULL) {
        cout << "Constructing coordinate system from file " << filename << endl;
        ctx.coords_file = new NcFile(filename, NcFile::read);
        string fs_dimensions;
        ctx.coords_file->getAtt("featurespace_dimensions").getValues(fs_dimensions);
        ctx.dim_names = vectors::from_string<string>(fs_dimensions);
        string fs_vars;
        ctx.coords_file->getAtt("featurespace_variables").getValues(fs_vars);
        vector<string> all_var_names = vectors::from_string<string>(fs_vars);
        for (size_t i = ctx.spatial_rank; i < all_var_names.size(); i++) {
            ctx.var_names.push_back(all_var_names[i]);
        }
#if WITH_VTK
        cout << "Setting VTK dimensions:" << ctx.params.vtk_dim_names << endl;
        VisitUtils<FS_TYPE>::update_vtk_dimension_mapping(ctx.dim_names, ctx.params.vtk_dim_names);
#endif
        // Construct coordinate system
        multimap<string, NcVar> vars = ctx.coords_file->getVars();
        multimap<string, NcVar>::iterator vi;
        cout << "Collating dimensions and dimension data ..." << endl;
        for (size_t di = 0; di < ctx.dim_names.size(); di++) {
            NcDim dim = ctx.coords_file->getDim(ctx.dim_names[di]);
            if (dim.isNull()) {
                cerr << "FATAL: dimension " << ctx.dim_names[di] << " does not exist " << endl;
                exit(EXIT_FAILURE);
            }
            NcVar var = ctx.coords_file->getVar(ctx.dim_names[di]);
            if (var.isNull()) {
                cerr << "FATAL: variable " << ctx.dim_names[di] << " does not exist " << endl;
                exit(EXIT_FAILURE);
            }
            ctx.dimensions.push_back(dim);
            ctx.dimension_vars.push_back(var);
        }
        cout << "Spatial range variables (dimensions): " << ctx.dim_names << endl;
        cout << "Value range variables: " << ctx.var_names << endl;
        ctx.coord_system = new CoordinateSystem<FS_TYPE>(ctx.dimensions, ctx.dimension_vars);
    }
}

template <typename T>
void readTrackingData(trackstats_context_t &ctx) {
    // run the directory's files
    fs::directory_iterator dir_iter(ctx.params.sourcepath);
    fs::directory_iterator end;
    while (dir_iter != end) {
        //            double time_reading = 0;
        //            double time_searching = 0;
        //            double time_inserting = 0;
        //            double time_pushing_clusters = 0;
        //            double time_constructing_clusters = 0;
        //            double time_deleting = 0;
        //            double time_constructing = 0;

        fs::path f = dir_iter->path();
        if (fs::is_regular_file(f) && boost::algorithm::ends_with(f.filename().generic_string(), "-clusters.nc")) {

            std::string path = f.generic_string();
            std::string filename = f.filename().generic_string();

            // Set up coordinate system, variable- and dimension names
            getFeaturespaceInfo(ctx, path);

            try {
                // start_timer();
                ctx.cluster_list = ClusterList<FS_TYPE>::read(path);
                int timeDifference = ctx.cluster_list->tracking_time_difference;
                // time_reading += stop_timer();

                cout << "Processing " << filename << " (" << ctx.cluster_list->size() << " clusters) ... ";

                if (ctx.spatial_rank == 0) {
                    ctx.spatial_rank = ctx.cluster_list->dimensions.size();
                } else if (ctx.spatial_rank != ctx.cluster_list->dimensions.size()) {
                    cerr << "FATAL:spatial range must remain identical across the track" << endl;
                    exit(EXIT_FAILURE);
                }

                size_t v_rank = ctx.cluster_list->feature_variables.size() - ctx.spatial_rank;
                if (ctx.value_rank == 0) {
                    ctx.value_rank = v_rank;
                } else if (v_rank != ctx.value_rank) {
                    cerr << "FATAL:value range must remain identical across the track" << endl;
                    exit(EXIT_FAILURE);
                }

                // Iterate over the clusters in the list we just read
                Cluster<FS_TYPE>::list::const_iterator ci;
                for (ci = ctx.cluster_list->clusters.begin(); ci != ctx.cluster_list->clusters.end(); ++ci) {
                    Cluster<FS_TYPE>::ptr cluster = (*ci);
                    m3D::id_t id = id = cluster->id;
                    Track<FS_TYPE>::ptr tm = NULL;
                    Track<FS_TYPE>::trackmap::const_iterator ti;

                    // start_timer();
                    ti = ctx.track_map.find(id);
                    // time_searching += stop_timer();

                    if (ti == ctx.track_map.end()) {
                        // new entry
                        tm = new Track<FS_TYPE>();
                        tm->id = id;
                        // start_timer();
                        ctx.track_map[id] = tm;
                        // time_inserting += stop_timer();
                    } else {
                        tm = ti->second;
                    }

                    boost::filesystem::path sf(ctx.cluster_list->source_file);
                    tm->sourcefiles.push_back(sf.filename().generic_string());

                    // Instead of the original cluster, use a TrackCluster, which
                    // has facilities of writing it's point list to disk and read
                    // it back on demand, saving memory.

                    // start_timer();
                    TrackCluster<FS_TYPE>::ptr tc = new TrackCluster<FS_TYPE>(cluster,timeDifference,ctx.need_points);
                    // time_constructing_clusters += stop_timer();

                    node_t node;
                    node.uuid = cluster->uuid;
                    node.id = cluster->id;
                    node.step = ctx.step;
                    node.size = cluster->size();

                    addGraphNode(ctx, node);

                    // Keep track to calculate average cluster size
                    ctx.average_cluster_size += cluster->size();
                    ctx.number_of_clusters++;

                    // start_timer();
                    tm->clusters.push_back(tc);
                    // time_pushing_clusters += stop_timer();

                    // Pre-empt calculations based on points
                    // and dispose of the data if possible

                    tc->geometrical_center();
                    vector<FS_TYPE> min, max, median;
                    cluster->variable_ranges(min, max, median);
                    ctx.cluster_min[cluster->id] = min;
                    ctx.cluster_max[cluster->id] = max;
                    ctx.cluster_median[cluster->id] = median;
                    // start_timer();
                    cluster->clear(true);
                    // time_deleting += stop_timer();
                }

                // start_timer();
                delete ctx.cluster_list;
                ctx.cluster_list = NULL;
                // time_deleting += stop_timer();

                // cout << "tracking map has " << track_map.size() << " tracks" << endl;
                cout << "done." << endl;

            } catch (netCDF::exceptions::NcException &e) {
                cerr << "ERROR:" << e.what() << endl;
            } catch (std::exception &e) {
                cerr << "ERROR:" << e.what() << endl;
            }
        }

        // cout << "Clusters processed:" << endl;
        // cout << "\t number:" << number_of_clusters << endl;
        // cout << "\t average size:" << ((double)average_cluster_size) / ((double)number_of_clusters) << endl;
        // cout << endl;
        //
        // cout << "Time spend:" << endl;
        // cout << "\t reading:" << time_reading << endl;
        // cout << "\t searching:" << time_searching << endl;
        // cout << "\t inserting:" << time_inserting << endl;
        // cout << "\t constructing clusters:" << time_constructing_clusters << endl;
        // cout << "\t pushing clusters:" << time_pushing_clusters << endl;
        // cout << "\t deleting:" << time_deleting << endl;

        dir_iter++;
        ctx.step++;
    }
    
    // Clear out all point data
    Track<FS_TYPE>::trackmap::iterator tmi;
    for (tmi = ctx.track_map.begin(); tmi != ctx.track_map.end(); ++tmi) {
        Track<FS_TYPE>::ptr track = tmi->second;
        size_t i = 0;
        std::list<Cluster<FS_TYPE>::ptr>::const_iterator ti;
        for (ti = track->clusters.begin(); ti != track->clusters.end(); ++ti) {
            TrackCluster<FS_TYPE> *c = (TrackCluster<FS_TYPE> *) * ti;
            c->clear(true);
        }
    }
}

#pragma mark -
#pragma mark Data collation

/**
 * Gets the speed histogram data from the currently
 * processed track.
 * 
 * @param ctx
 */
void getSpeedHistogramData(trackstats_context_t &ctx) {
    FS_TYPE dS = boost::numeric_cast<FS_TYPE>(1000.0) * vector_norm(ctx.cluster->displacement);
    if (dS == 0) return;
    FS_TYPE dT = boost::numeric_cast<FS_TYPE>(ctx.cluster->tracking_time_difference());
    FS_TYPE speed = dS / dT;
    bool exceeded_max_class = false;
    add_value_to_histogram<float>(ctx.params.speed_histogram_bins,
            ctx.speed_histogram,
            speed,
            exceeded_max_class);
    if (exceeded_max_class) {
        cout << "Cluster uuid:" << ctx.cluster->uuid
             << " id:" << ctx.cluster->id
             << " (speed " << speed << " m/s)"
             << " exceeded max speed "
             << ctx.params.speed_histogram_bins.back() << " m/s"
             << endl;
    }
    map<float, size_t>::iterator si = ctx.speeds.find(speed);
    if (si == ctx.speeds.end()) {
        ctx.speeds[speed] = 1;
    } else {
        ctx.speeds[speed] = (si->second + 1);
    }
}

/**
 * Gets the direction histogram data from the currently
 * processed track.
 * 
 * TODO: This code is using the standard nationwide composite here
 * and so is dependant on libradolan and fixed to the RADOLAN
 * coordinate system which is BAD. Find a more generic way!!
 * 
 * @param ctx
 */
void getDirectionHistogramData(trackstats_context_t &ctx) {
    RDCoordinateSystem *rcs = new RDCoordinateSystem(RD_RX);

    // calculate speed in m/s. All vectors 2D at this time
    vector<FS_TYPE> dP = reinterpret_cast<Cluster<FS_TYPE>::ptr>(ctx.cluster)->displacement;
    if (vector_norm(dP) == 0) return;

    // Assuming --vtk-dimensions=x,y
    RDCartesianPoint p;
    p.x = ctx.cluster->displacement.at(1);
    p.y = ctx.cluster->displacement.at(0);

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
    if (direction > 360.0) {
        direction -= 360.0;
    }

    bool exceeded_max_class = false;
    add_value_to_histogram<float>(ctx.params.direction_histogram_bins,
            ctx.direction_histogram,
            direction,
            exceeded_max_class);
    map<float, size_t>::iterator si = ctx.directions.find(direction);
    if (si == ctx.directions.end()) {
        ctx.directions[direction] = 1;
    } else {
        ctx.directions[direction] = (si->second + 1);
    }
}

/**
 * Gets the size histogram data from the currently
 * processed track.
 * 
 * @param ctx
 */
void getSizeHistogramData(trackstats_context_t &ctx) {
    size_t cluster_size = ctx.cluster->size();
    bool exceeded_max_class = false;
    add_value_to_histogram(ctx.params.cluster_histogram_bins, 
            ctx.cluster_histogram, 
            cluster_size, 
            exceeded_max_class);
    if (exceeded_max_class) {
        cout << "Cluster uuid:" << ctx.cluster->uuid
                << " id:" << ctx.cluster->id
                << " (size " << cluster_size << ")"
                << " exceeded max cluster size"
                << ctx.params.cluster_histogram_bins.back() << endl;
    }
    map<size_t, size_t>::iterator csfi = ctx.cluster_sizes.find(cluster_size);
    if (csfi == ctx.cluster_sizes.end()) {
        ctx.cluster_sizes[cluster_size] = 1;
    } else {
        ctx.cluster_sizes[cluster_size] = (csfi->second + 1);
    }
}

/**
 * Iterate over the points of the cluster and add points to the 
 * array index. If the point exists, add it's values to the existing 
 * one. If not, add a new point. In the end the index contains all 
 * points that all clusters in the track occupied at any time
 *
 * @param ctx
 */
void addToCumulativeStats(trackstats_context_t &ctx) {
    Point<FS_TYPE>::list::iterator pi;
    for (pi = ctx.cluster->get_points().begin(); 
            pi != ctx.cluster->get_points().end(); ++pi) 
    {
        Point<FS_TYPE>::ptr p = *pi;
        Point<FS_TYPE>::ptr indexed = ctx.track_index->get(p->gridpoint);
        if (indexed == NULL) {
            // this makes a copy of the point in the index
            ctx.track_index->set(p->gridpoint, p);
            // get the copied point
            indexed = ctx.track_index->get(p->gridpoint);
        }
        // only add up the value range
        for (size_t k = 0; k < ctx.value_rank; k++) {
            indexed->values[ctx.spatial_rank + k] += p->values[ctx.spatial_rank + k];
        }
        
        ctx.points_processed++;
    }
}

/**
 * Gets the size histogram data from the currently
 * processed track.
 * 
 * @param ctx
 */
void getLengthHistogramData(trackstats_context_t &ctx) {
    size_t track_length = ctx.track->size();
    bool exceeded_max_class = false;
    add_value_to_histogram(ctx.params.length_histogram_bins, 
            ctx.length_histogram, 
            track_length, 
            exceeded_max_class);
    if (exceeded_max_class) {
        cout << "Track id:" << ctx.track->id
        << " (track length " << track_length << " m/s)"
        << " exceeded max speed "
        << ctx.params.speed_histogram_bins.back() << " m/s"
        << endl;
    }

    map<size_t, size_t>::iterator tlfi = ctx.track_lengths.find(track_length);
    if (tlfi == ctx.track_lengths.end())
        ctx.track_lengths[track_length] = 1;
    else
        ctx.track_lengths[track_length] = (tlfi->second + 1);
}

void getCumulatedStatsData(trackstats_context_t& ctx) {
    size_t cumulatedSize = ctx.track_index->count();
    bool exceeded_max_class = false;
    add_value_to_histogram(ctx.params.size_histogram_bins, 
            ctx.size_histogram, 
            cumulatedSize, 
            exceeded_max_class);

    if (exceeded_max_class) {
        cout << "Track id:" << ctx.track->id
        << " (cumulative size " << cumulatedSize
        << " exceeded max cumulated size class: "
        << ctx.params.size_histogram_bins.back()
        << endl;
    }

    map<size_t, size_t>::iterator tlfi = ctx.track_sizes.find(cumulatedSize);
    if (tlfi == ctx.track_sizes.end()) {
        ctx.track_sizes[cumulatedSize] = 1;
    } else {
        ctx.track_sizes[cumulatedSize] = (tlfi->second + 1);
    }
}

#pragma mark -
#pragma mark Writing outputs

void writeTrackDictionary(const trackstats_context_t &ctx) {
    ofstream dict("track-dictionary.json");

    dict << "{" << endl;
    dict << "  \"spatial_range\":" << to_json(ctx.dim_names) << "," << endl;
    dict << "  \"value_range\":" << to_json(ctx.var_names) << "," << endl;
    dict << "  \"number_of_tracks\":" << ctx.track_map.size() << "," << endl;

    dict << "  " << "\"tracks\":[" << endl;

    size_t j = 0;
    Track<FS_TYPE>::trackmap::const_iterator tmi;
    for (tmi = ctx.track_map.begin(); tmi != ctx.track_map.end(); ++tmi) {
        Track<FS_TYPE>::ptr track = tmi->second;

        // Exclude degenerates if that's in the cards
        if (track->clusters.size() == 1 && ctx.params.exclude_degenerates) {
            continue;
        }

        if (j > 0 && j < (ctx.track_map.size()-1)) {
            dict << "," << endl;
        }

        dict << "    {" << endl;
        dict << "      \"id\":" << track->id << "," << endl;
        dict << "      \"length\":" << track->clusters.size() << "," << endl;
        dict << "      \"clusters\":[" << endl;

        track->min.resize(ctx.value_rank, std::numeric_limits<FS_TYPE>::max());
        track->max.resize(ctx.value_rank, std::numeric_limits<FS_TYPE>::min());

        size_t i = 0;
        std::list<Cluster<FS_TYPE>::ptr>::const_iterator ti;
        for (ti = track->clusters.begin(); ti != track->clusters.end(); ++ti) {
            TrackCluster<FS_TYPE> *c = (TrackCluster<FS_TYPE> *) (*ti);
            vector<FS_TYPE> min, max, median;
            min = ctx.cluster_min.at(c->id);
            max = ctx.cluster_max.at(c->id);
            median = ctx.cluster_median.at(c->id);

            dict << "        {" << endl;
            dict << "          \"size\":" << c->size() << "," << endl;
            dict << "          \"uuid\":" << c->uuid << "," << endl;
            dict << "          \"sourcefile\":\"" << track->sourcefiles[i] << "\"," << endl;
            dict << "          \"geometrical_center\":" << to_json(c->geometrical_center()) << "," << endl;
            dict << "          \"mode\":" << to_json(c->mode) << "," << endl;
            dict << "          \"min\":" << to_json(min) << "," << endl;
            dict << "          \"max\":" << to_json(max) << "," << endl;
            dict << "          \"median\":" << to_json(median) << "," << endl;
            dict << "          \"has_margin_points\":" << (c->has_margin_points() ? "true" : "false") << endl;
            dict << "        }";

            // min/max
            for (size_t k = 0; k < ctx.value_rank; k++) {
                if (min[k] < track->min[k]) track->min[k] = min[k];
                if (max[k] > track->max[k]) track->max[k] = max[k];
            }

            if (i < (track->clusters.size() - 1))
                dict << "," << endl;

            // Dispose of the points again to free up memory.
            c->clear(true);

            i++;
        }
        dict << endl << "      ]," << endl;

        // Limits / Median
        dict << "      \"min\":" << to_json(track->min) << "," << endl;
        dict << "      \"max\":" << to_json(track->max) << endl;
        dict << "    }";

        j++;
    }

    // end tracks
    dict << "  ],";

    // begin tree data
    dict << "  " << "\"tree\":{" << endl;

    // begin nodes
    dict << "    " << "\"nodes\":[" << endl;
    for (size_t ni = 0; ni < ctx.nodes.size(); ni++) {
        node_t node = ctx.nodes[ni];
        dict << "      {" << endl;
        dict << "        \"uuid\":" << node.uuid << "," << endl;
        dict << "        \"id\":" << node.id << "," << endl;
        dict << "        \"size\":" << node.size << "," << endl;
        dict << "        \"step\":" << node.step << endl;
        dict << "      }";
        if (ni < (ctx.nodes.size() - 1)) {
            dict << "," << endl;
        }
    }
    // end nodes
    dict << "    " << "]," << endl;

    // begin links
    dict << "    " << "\"links\":[" << endl;
    for (size_t ni = 0; ni < ctx.links.size(); ni++) {
        link_t link = ctx.links[ni];
        dict << "      {" << endl;
        dict << "        \"source\":" << link.source << "," << endl;
        dict << "        \"target\":" << link.target << "," << endl;
        dict << "        \"type\":" << link.type << "," << endl;
        dict << "        \"id\":" << link.id << endl;
        dict << "      }";
        if (ni < (ctx.links.size() - 1)) {
            dict << "," << endl;
        }
    }
    // end links
    dict << "    " << "]" << endl;

    // end tree data
    dict << "  " << "}" << endl;

    // end dictionary
    dict << "}" << endl;
    dict.close();
}

void writeStats(trackstats_context_t &ctx) {
    string file_fn = ctx.params.basename + "_trackstats.txt";
    ofstream file(file_fn.c_str());

    // number of tracks
    file << "Tracking report for basename " + ctx.params.basename << endl;
    cout << "Tracking report for basename " + ctx.params.basename << endl;

    file << "(degenerates are excluded: " << (ctx.params.exclude_degenerates ? "no" : "yes") << ")" << endl;
    cout << "(degenerates are excluded: " << (ctx.params.exclude_degenerates ? "no" : "yes") << ")" << endl;

    file << endl;
    cout << endl;

    file << "Overall number of tracks: " << ctx.track_map.size() << endl;
    cout << "Overall number of tracks: " << ctx.track_map.size() << endl;

    file << "Number of degenerate tracks: " << ctx.number_of_degenerates << endl;
    cout << "Number of degenerate tracks: " << ctx.number_of_degenerates << endl;

    file << endl;
    cout << endl;
    
    // length of tracks (in time steps)

    if (ctx.params.create_length_stats) {
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Track length distribution" << endl;
        cout << "Track length distribution" << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<size_t>(ctx.track_lengths, true);
        file << "(Average track length = " << avg << ")" << endl;
        cout << "(Average track length = " << avg << ")" << endl;

        map<size_t, size_t>::reverse_iterator rend = ctx.track_lengths.rbegin();
        file << "(Maximum track length = " << rend->first << ")" << endl;
        cout << "(Maximum track length = " << rend->first << ")" << endl;

        file << "length,number" << endl;
        cout << "length,number" << endl;

        map<size_t, size_t>::iterator si;

        for (si = ctx.track_lengths.begin(); 
                si != ctx.track_lengths.end(); si++) {
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

        print_histogram(ctx.params.length_histogram_bins, ctx.length_histogram, file);

        if (ctx.params.write_gnuplot_files) {
            write_histogram("lengths-hist.txt", "number of tracks", "track length", 
                    ctx.params.length_histogram_bins, 
                    ctx.length_histogram);
            
            write_values<size_t>("lengths.txt", "number of tracks", "track length", 
                    ctx.track_lengths);
        }
    }

    if (ctx.params.create_speed_stats) {
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Speed distribution" << endl;
        cout << "Speed distribution" << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<float>(ctx.speeds, true);
        file << "(Average speed = " << avg << ")" << endl;
        cout << "(Average speed = " << avg << ")" << endl;

        map<float, size_t>::reverse_iterator srend = ctx.speeds.rbegin();
        file << "(Maximum speed = " << srend->first << ")" << endl;
        cout << "(Maximum speed = " << srend->first << ")" << endl;

        file << "speed,number" << endl;
        cout << "speed,number" << endl;

        map<float, size_t>::iterator spi;

        for (spi = ctx.speeds.begin(); spi != ctx.speeds.end(); spi++) {
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

        print_histogram<float>(ctx.params.speed_histogram_bins, 
                ctx.speed_histogram, file);

        if (ctx.params.write_gnuplot_files) {
            write_histogram("speeds-hist.txt", "number of clusters", "speed [m/s]", 
                    ctx.params.size_histogram_bins, 
                    ctx.speed_histogram);
            write_values<float>("speeds.txt", "number of clusters", "speed in [m/s]", 
                    ctx.speeds);
        }
    }

    if (ctx.params.create_direction_stats) {
        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Direction distribution" << endl;
        cout << "Direction distribution" << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<float>(ctx.directions, true);
        file << "(Average direction = " << avg << ")" << endl;
        cout << "(Average direction = " << avg << ")" << endl;

        file << "direction,number" << endl;
        cout << "direction,number" << endl;

        map<float, size_t>::iterator spi;

        for (spi = ctx.directions.begin(); 
                spi != ctx.directions.end(); spi++) {
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

        print_histogram<float>(ctx.params.direction_histogram_bins, 
                ctx.direction_histogram, file);

        if (ctx.params.write_gnuplot_files) {
            write_histogram("directions-hist.txt", "number of clusters", "direction in [deg]", 
                    ctx.params.direction_histogram_bins, 
                    ctx.direction_histogram);
            write_values<float>("directions.txt", "number of clusters", "direction in [deg]", 
                    ctx.directions);
        }
    }

    if (ctx.params.create_cluster_stats) {
        file << endl;
        cout << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Cluster size distribution" << endl;
        cout << "Cluster size distribution" << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<size_t>(ctx.cluster_sizes, true);
        file << "(Average cluster size = " << avg << ")" << endl;
        cout << "(Average cluster size = " << avg << ")" << endl;

        map<size_t, size_t>::reverse_iterator rend = ctx.cluster_sizes.rbegin();
        file << "(Maximum cluster size = " << rend->first << ")" << endl;
        cout << "(Maximum cluster size = " << rend->first << ")" << endl;

        file << "size,number" << endl;
        cout << "size,number" << endl;

        map<size_t, size_t>::iterator si;

        for (si = ctx.cluster_sizes.begin(); si != ctx.cluster_sizes.end(); si++) {
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

        print_histogram(ctx.params.cluster_histogram_bins, 
                ctx.cluster_histogram, file);

        if (ctx.params.write_gnuplot_files) {
            write_histogram("sizes-hist.txt", "number of clusters", "size [#gridpoints]", 
                    ctx.params.cluster_histogram_bins, 
                    ctx.cluster_histogram);
            write_values("sizes.txt", "number of clusters", "size [#gridpoints]", 
                    ctx.cluster_sizes);
        }
    }

    if (ctx.params.create_cumulated_size_stats) {
        // size of tracks (in terms of cumulative pixels)

        file << endl;
        cout << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        file << "Track size distribution" << endl;
        cout << "Track size distribution" << endl;

        file << "------------------------------------------------" << endl;
        cout << "------------------------------------------------" << endl;

        double avg = average<size_t>(ctx.track_sizes, false);
        file << "(Average track size = " << avg << ")" << endl;
        cout << "(Average track size = " << avg << ")" << endl;

        map<size_t, size_t>::reverse_iterator rend = ctx.track_sizes.rbegin();
        file << "(Maximum track size = " << rend->first << ")" << endl;
        cout << "(Maximum track size = " << rend->first << ")" << endl;

        file << "size,number" << endl;
        cout << "size,number" << endl;

        map<size_t, size_t>::iterator si;

        for (si = ctx.track_sizes.begin(); si != ctx.track_sizes.end(); si++) {
            file << si->first << "," << si->second << endl;
            cout << si->first << "," << si->second << endl;
        }

        file << endl;
        cout << endl;

        file << "Size Histogram:" << endl;
        cout << "Size Histogram:" << endl;

        file << "max size,number" << endl;
        cout << "max size,number" << endl;

        print_histogram(ctx.params.size_histogram_bins, ctx.size_histogram, file);

        if (ctx.params.write_gnuplot_files) {
            write_histogram("cumulative-sizes-hist.txt", "number of tracks", "size [#gridpoints]", 
                    ctx.params.size_histogram_bins, ctx.size_histogram);
            write_values<size_t>("cumulative-sizes.txt", "number of tracks", "size [#gridpoints]", 
                    ctx.track_sizes);
        }
    }
}

# pragma mark -
# pragma mark Processing control

/**
 * Processes the entire track map.
 * 
 * @param ctx
 */
void processTracks(trackstats_context_t &ctx) {
    
    cout << endl << "Keying up tracks: " << endl;

    // Iterate over the collated tracks
    Track<FS_TYPE>::trackmap::iterator tmi;
    for (tmi = ctx.track_map.begin(); tmi != ctx.track_map.end(); ++tmi) {

        ctx.track = tmi->second;
        
        // Handle degenerate tracks?
        if (ctx.track->clusters.size() == 1) {
            ctx.number_of_degenerates++;
            if (ctx.params.exclude_degenerates) continue;
        }

        cout << "Processing track #" << tmi->first 
             <<  " (" << ctx.track->clusters.size() << " clusters)" 
             << endl;

        if (ctx.params.create_length_stats) {
            getLengthHistogramData(ctx);
        }

        // Create an array index for this track. Start with empty index.
        if (ctx.params.create_cumulated_size_stats) {
            ctx.track_index = new ArrayIndex<FS_TYPE>(ctx.coord_system->get_dimension_sizes(), false);
        }

        // Process the track's clusters
        std::list<Cluster<FS_TYPE>::ptr>::iterator ti;
        for (ti = ctx.track->clusters.begin(); 
                ti != ctx.track->clusters.end(); ++ti) 
        {
            ctx.cluster = (TrackCluster<FS_TYPE> *) *ti;

            if (ctx.params.create_cluster_stats) {
                getSizeHistogramData(ctx);
            }

            if (ctx.params.create_speed_stats) {
                getSpeedHistogramData(ctx);
            }

            // Directional Statistics
            if (ctx.params.create_direction_stats) {
                getDirectionHistogramData(ctx);
            }

            if (ctx.params.create_cumulated_size_stats) {
                addToCumulativeStats(ctx);
            }
            
            // Purge points in memory
            ctx.previous_cluster = ctx.cluster;
            ctx.cluster->clear(true);
        }

#if WITH_VTK
        if (ctx.params.write_cumulated_tracks_as_vtk) {
            Point<FS_TYPE>::list cumulatedList;
            ctx.track_index->replace_points(cumulatedList);
            string vtk_path = ctx.params.basename + "_cumulated_track_" 
                    + boost::lexical_cast<string>(tmi->first) + ".vtk";
            VisitUtils<FS_TYPE>::write_pointlist_all_vars_vtk(vtk_path, 
                    &cumulatedList, vector<string>());
        }
#endif
        // Delete the index to save memory
        if (ctx.params.create_cumulated_size_stats) {
            // Levy the stats from the cumulated data
            getCumulatedStatsData(ctx);
            // Delete the index
            ctx.track_index->clear(true);
            delete ctx.track_index;
            ctx.track_index = NULL;
            
            cout << "  (processed " << ctx.points_processed << " points)" << endl;
        }
        delete ctx.track;
    }
    cout << "done." << endl;
}

int main(int argc, char** argv) {
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
            ("sourcepath,s", program_options::value<string>()->default_value("."), "Directory containing cluster files.")
            ("exclude-degenerates,n", program_options::value<bool>()->default_value(true), "Exclude results of tracks of length one")
            ("create-length-statistics,1", "Create a statistic of track lengths.")
            ("length-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(length_hist_default), "List of track-length values for histogram bins")
            ("create-speed-statistics,2", "Evaluate speeds of clusters in tracks, based on geometric center point")
            ("speed-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(speed_hist_default), "Speed histogram. Values in [m/s]")
            ("create-direction-statistics,3", "Evaluate directions of clusters in tracks, based on geometric center point")
            ("direction-histogram-classes", program_options::value<fvec_t>()->multitoken()->default_value(direction_hist_default), "Direction histogram. Values in [deg]. Use with radolan grid only!!")
            ("create-cluster-statistics,4", "Evaluate each cluster in each track in terms of size.")
            ("cluster-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(cluster_hist_default), "List of cluster size values for histogram bins")
            ("create-cumulated-tracking-stats,5", "Obtain a number of tracking stats")
            ("create-cumulated-size-statistics,6", "Evaluate each track in terms of cumulative size. Warning: this process takes a lot of memory.")
            ("size-histogram-classes", program_options::value<bin_t>()->multitoken()->default_value(size_hist_default), "List of cumulated track size values for histogram bins")
            ("write-track-dictionary,t", "Write out a dictionary listing tracks with number of clusters etc.")
#if WITH_VTK
            ("write-center-tracks-as-vtk,e", "Write tracks out as .vtk files")
            ("write-cumulated-tracks-as-vtk,m", "Write cumulated tracks out as .vtk files. Only has effect if --create-cumulated-size-statistics is used")
            ("vtk-dimensions", program_options::value<string>(), "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)")
#endif    
            ("write-gnuplot-files,g", "write individual files for the statistics fit for use with gnuplot")
            ;

    program_options::variables_map vm;

    try {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    } catch (std::exception &e) {
        cerr << "FATAL:Error parsing command line: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("help") == 1 || argc < 2) {
        cout << desc << "\n";
        return EXIT_SUCCESS;
    }

    if (vm.count("version") != 0) {
        cout << m3D::VERSION << endl;
        return EXIT_SUCCESS;
    }

    // Evaluate user input
    parameter_t params;
    try {
        parse_commmandline(vm, params);
    } catch (const std::exception &e) {
        cerr << "FATAL:" << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // initialize the run context data
    trackstats_context_t ctx = initialiseContext(params);

    // This file should not be closed until the code has
    // run through, or netCDF will upchuck exceptions when
    // accessing coordinate system
    if (!fs::is_directory(ctx.params.sourcepath)) {
        cerr << "Argument --sourcepath does not point to a directory "
                "(--sourcepath=" << ctx.params.sourcepath << ")" 
                << endl;
        return EXIT_FAILURE;
    }

    readTrackingData<FS_TYPE>(ctx);

#if WITH_VTK
    if (ctx.params.write_center_tracks_as_vtk) {
        // Write center tracks
        VisitUtils<FS_TYPE>::write_center_tracks_vtk(ctx.track_map,
                ctx.params.basename,
                ctx.spatial_rank,
                ctx.params.exclude_degenerates);
    }
#endif        
    // dictionary?

    if (ctx.params.write_track_dictionary) {
        writeTrackDictionary(ctx);
    }

    processTracks(ctx);
    writeStats(ctx);

    // Clean up
    delete ctx.coords_file;

    // Done.
    return EXIT_SUCCESS;
};
