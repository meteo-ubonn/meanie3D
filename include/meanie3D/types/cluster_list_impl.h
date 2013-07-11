#ifndef _M3D_ClusterList_Impl_H_
#define _M3D_ClusterList_Impl_H_

#include <algorithm>
#include <sstream>
#include <netcdf>
#include <vector>
#include <set>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <stdlib.h>

namespace m3D {

	using namespace std;
    using namespace netCDF;
    using namespace cfa::utils::console;
    using namespace cfa::utils::timer;
    using namespace cfa::utils::vectors;
    using namespace cfa::utils::visit;
	using cfa::meanshift::Point;
	using cfa::meanshift::FeatureSpace;
	using cfa::meanshift::PointIndex;
	using cfa::meanshift::KNNSearchParams;
	using cfa::meanshift::RangeSearchParams;
	using cfa::meanshift::SearchParameters;

#pragma mark -
#pragma mark Macros
    
	#define sup( v1,v2 ) (v1 > v2 ? v1:v2)
	#define inf( v1,v2 ) (v1 < v2 ? v1:v2)
    
#pragma mark -
#pragma mark Accessing the list

    template <typename T>
    size_t
    ClusterList<T>::size()
    {
        return clusters.size();
    }
    
    template <typename T>
    typename Cluster<T>::ptr
    ClusterList<T>::operator[] (size_t index)
    {
        return clusters[index];
    }

    
#pragma mark -
#pragma mark Adding / Removing points
    
    template <typename T>
    void
    ClusterList<T>::apply_size_threshold( unsigned int min_cluster_size, const bool& show_progress )
    {
        boost::progress_display *progress = NULL;
        
        if ( show_progress )
        {
            cout << endl << "Applying size threshold of " << min_cluster_size << " ... " << endl;
            
            start_timer();
            
            progress = new boost::progress_display(clusters.size());
        }
        
        size_t axe_count = 0;
        
        if ( min_cluster_size > 1 )
        {
            typename Cluster<T>::list::iterator it;
            
            for ( it = clusters.begin(); it != clusters.end(); )
            {
                progress->operator++();
                
                typename Cluster<T>::ptr sc = *it;
                
                if ( sc->points.size() < min_cluster_size )
                {
                    it = clusters.erase( it );
                    
                    axe_count++;
                }
                else
                {
                    it++;
                }
            }
        }
        
        if ( show_progress )
        {
            cout << endl << "done. (Removed " << axe_count << " objects in " << stop_timer() << "s)" << endl;
            
            delete progress;
        }
    }
    
    template <typename T>
    void 
    ClusterList<T>::write( const std::string& path )
    {
        try
        {
            // use NetCDF C - API (instead of C++) because the C++ API has no
            // support for variable length arrays :(
            
            bool fileExists = this->ncFile != NULL;
            
            NcFile *file = NULL;
            
            try
            {
                if ( fileExists )
                {
                    file = this->ncFile;
                }
                else
                {
                    file = new NcFile( path, NcFile::replace );
                }
            }
            catch ( const netCDF::exceptions::NcException &e )
            {
                cerr << "Exception opening file " << path << " for writing : " << e.what() << endl;
                exit(-1);
            }
            
            // Create dimensions
            
            NcDim dim,spatial_dim;
            
            if ( fileExists )
            {
                dim = file->getDim("featurespace_dim");
                
                spatial_dim = file->getDim("spatial_dim");
            }
            else
            {
                dim = file->addDim("featurespace_dim", (int) this->feature_variables.size() );
            
                spatial_dim = file->addDim("spatial_dim", (int) this->dimensions.size() );
                
                // write featurespace_dimensions attribute
                
                vector<string> fs_dims;
                
                for (size_t di=0; di<this->dimensions.size(); di++)
                {
                    fs_dims.push_back(this->dimensions[di].getName());
                }
                
                file->putAtt("featurespace_dimensions", to_string(fs_dims) );
                
                // copy dimensions
                
                for (size_t di=0; di<this->dimensions.size(); di++)
                {
                    NcDim d = this->dimensions[di];
                    
                    file->addDim(d.getName(), d.getSize());
                }
                
                // Create dummy variables, attributes and other meta-info
                
                file->putAtt( "num_clusters", ncInt, (int) clusters.size() );
                
                file->putAtt( "source", this->source_file );

                // compile a list of the feature space variables used,
                // including the spatial dimensions
                
                vector<string> variable_names;
                
                // feature variables
                
                for ( size_t i=0; i < this->feature_variables.size(); i++ )
                {
                    NcVar var = this->feature_variables[i];

                    // append to list
                    
                    variable_names.push_back(var.getName());
                    
                    NcDim *dim = NULL;

                    // Copy data (in case of dimension variables)
                    // exploiting once more the fact, that dimension variables
                    // have by convention the same name as the dimension
                    
                    for (size_t di=0; di<this->dimensions.size() && dim==NULL; di++)
                    {
                        NcDim d = this->dimensions[di];
                        dim = (d.getName() == var.getName()) ? &d : NULL;
                    }
                    
                    // create a dummy variable
                    
                    NcVar dummyVar;
                    
                    if (dim != NULL)
                    {
                        dummyVar = file->addVar( var.getName(), var.getType(), *dim);
                        
                        T *data = (T*)malloc(sizeof(T) * dim->getSize());
                        
                        var.getVar(data);
                        
                        dummyVar.putVar(data);
                        
                        delete data;
                    }
                    else
                    {
                        dummyVar = file->addVar( var.getName(), var.getType(), spatial_dim );
                    }
                    
                    // Copy attributes
                    
                    map< string, NcVarAtt > attributes = var.getAtts();
                    
                    map< string, NcVarAtt >::iterator at;
                    
                    for ( at = attributes.begin(); at != attributes.end(); at++ )
                    {
                        NcVarAtt a = at->second;
                        
                        size_t size = a.getAttLength();
                        
                        void *data = (void *)malloc( size );
                        
                        a.getValues( data );
                        
                        NcVarAtt copy = dummyVar.putAtt( a.getName(), a.getType(), size, data );
                        
                        free(data);
                    }
                    
                }
                
                file->putAtt("featurespace_variables", to_string(variable_names));
            }
            
            // Add tracking meta-info
            
            if ( this->tracking_performed )
            {
                file->putAtt( "tracking_performed", "yes" );
                file->putAtt( "tracked_ids", to_string( this->tracked_ids ) );
                file->putAtt( "new_ids", to_string( this->new_ids ) );
                file->putAtt( "dropped_ids", to_string( this->dropped_ids ) );
            }
            
            // Add cluster dimensions and variables
            
            for ( size_t ci = 0; ci < clusters.size(); ci++ )
            {
                // Create a dimension
                
                stringstream dim_name(stringstream::in | stringstream::out);
                
                dim_name << "cluster_dim_" << ci;
                
                NcDim cluster_dim;
                
                if ( fileExists )
                {
                    cluster_dim = file->getDim( dim_name.str() );
                }
                else
                {
                    cluster_dim = file->addDim( dim_name.str(), clusters[ci]->points.size() );
                }
                
                // Create variable
                
                stringstream var_name(stringstream::in | stringstream::out);
                
                var_name << "cluster_" << ci;
                
                vector<NcDim> dims(2);
                dims[0] = cluster_dim;
                dims[1] = dim;
                
                NcVar var;
                
                if ( fileExists )
                {
                    var = file->getVar( var_name.str() );
                }
                else
                {
                    var = file->addVar( var_name.str(), ncDouble, dims );
                }
                
                // size
                
                var.putAtt( "size", ncInt, (int) clusters[ci]->points.size() );
                
                // id
                
                unsigned long long cid = (unsigned long long) clusters[ci]->id;
                
                var.putAtt( "id", ncInt64, cid );
                
                // mode

                string mode = to_string( clusters[ci]->mode );

                var.putAtt( "mode", mode );

                // Write the clusters away point by point
                
                // index and counter
                vector<size_t> index(2,0);
                vector<size_t> count(2,0);
                count[0] = 1;
                count[1] = dim.getSize();

                // iterate over points
                for ( size_t pi = 0; pi < clusters[ci]->points.size(); pi++ )
                {
                    Point<T> *p = clusters[ci]->points[pi];
                    
                    double data[ dim.getSize() ];
                    
                    for ( size_t di = 0; di < dim.getSize(); di++ )
                    {
                        data[di] = (double)p->values[di];
                    }
                    index[0] = pi;
                    var.putVar(index, count, &data[0] );
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception while writing cluster file: " << e.what() << endl;
        }
    }

    template <typename T> 
    typename ClusterList<T>::ptr
    ClusterList<T>::read(const std::string& path)
    {
        // meta-info
        
        vector<NcVar>                   feature_variables;
        vector<NcDim>                   dimensions;
        string                          source_file;
        typename Cluster<T>::list       list;
        NcFile                          *file = NULL;

        // TODO: let this exception go up
        try
        {
            file = new NcFile( path, NcFile::write );
        }
        catch ( const netCDF::exceptions::NcException &e )
        {
            cerr << "Could not open file " << path << " : " << e.what() << endl;
            exit( -1 );
        }
        
        try
        {
            // Read the dimensions

            NcDim fs_dim = file->getDim( "featurespace_dim" );
            
            // Read the dimensions attribute
            
            string dim_str;
            
            file->getAtt("featurespace_dimensions").getValues(dim_str);
            
            // Fill the dimensions vector from that
            
            vector<string> fs_dim_names = from_string<string>(dim_str);
            
            for (size_t di=0; di < fs_dim_names.size(); di++)
            {
                NcDim d = file->getDim(fs_dim_names[di]);
                
                dimensions.push_back(d);
            }

            // Read global attributes
            
            file->getAtt("source").getValues( source_file );
            
            int number_of_clusters;
            
            file->getAtt("num_clusters").getValues( &number_of_clusters );
            
            // Read the feature-variables
            
            multimap<string,NcVar> vars = file->getVars();
            
            typename multimap<string,NcVar>::iterator vi;
            
            for ( vi = vars.begin(); vi != vars.end(); vi++ )
            {
                if (! boost::starts_with( vi->first, "cluster_"))
                {
                    feature_variables.push_back( vi->second );
                }
            }
            
            // Read clusters one by one
            
            for ( size_t i = 0; i < number_of_clusters; i++ )
            {
                // cluster dimension
                
                stringstream dim_name(stringstream::in | stringstream::out);
                
                dim_name << "cluster_dim_" << i;
                
                NcDim cluster_dim = file->getDim( dim_name.str().c_str() );
                
                // Read the variable
                
                stringstream var_name(stringstream::in | stringstream::out);
                
                var_name << "cluster_" << i;
                
                NcVar var = file->getVar( var_name.str().c_str() );
                
                // Decode mode
                
                typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
                
                boost::char_separator<char> sep(",");
                
                std::string mode_str;
                
                var.getAtt("mode").getValues(mode_str);
                
                boost::replace_all(mode_str, "(", "" );
                boost::replace_all(mode_str, ")", "" );
                
                tokenizer tokens( mode_str, sep );
                
                vector<T> mode;
                
                for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter )
                {
                    string token = *tok_iter;

                    mode.push_back( (T)atof( token.c_str() ) );
                }
                
                // Create a cluster object
                
                typename Cluster<T>::ptr cluster = new Cluster<T>( mode, dimensions.size() );
                
                // read cluster id
                
                var.getAtt("id").getValues( &(cluster->id) );
                
                // iterate over the data

                // Read the points, one by one
                
                vector<size_t> index(2,0);
                vector<size_t> count(2,0);
                count[0] = 1;
                count[1] = fs_dim.getSize();
                
                for ( size_t point_index = 0; point_index < cluster_dim.getSize(); point_index++ )
                {
                    // Allocate data
                    
                    T data[fs_dim.getSize()];
                    
                    // set index up and read
                    
                    index[0] = point_index;
                    
                    var.getVar(index,count,&data[0]);
                    
                    // copy data over to vectors
                    
                    vector<T> coordinate(dimensions.size(),0);
                    
                    vector<T> values(fs_dim.getSize(),0);
                    
                    for ( size_t i=0; i<fs_dim.getSize(); i++)
                    {
                        values[i] = data[i];
                        
                        if ( i < dimensions.size() )
                        {
                            coordinate[i] = values[i];
                        }
                    }
                    
                    // Create a point and add it to the cluster
                    
                    Point<T> *p = PointFactory<T>::get_instance()->create( coordinate, values );
                    
                    cluster->add_point(p);
                }
                
                list.push_back( cluster );
            }
        }
        catch (const std::exception &e)
        {
            cerr << e.what() << endl;
        }
        
        ClusterList<T>::ptr cl = new ClusterList<T>( list, dimensions, feature_variables, source_file );
        
        cl->ncFile = file;
        
        return cl;
    }
    
    template <typename T>
    void
    ClusterList<T>::print()
    {
        for ( size_t ci = 0; ci < clusters.size(); ci++ )
        {
            typename Cluster<T>::ptr c = clusters[ci];
            
            cout << "Cluster #" << ci << " at " << c->mode << " (" << c->points.size() << " points.)" << endl;
        }

    }
    
#pragma mark -
#pragma mark Clustering by Graph Theory
     
    template <typename T>
    void
    ClusterList<T>::find_neighbours(FeatureSpace<T> *fs,
                                    ArrayIndex<T> &arrayIndex,
                                    typename Point<T>::list &list,
                                    size_t dimensionIndex,
                                    typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = fs->coordinate_system->dimensions()[dimensionIndex];
        
        // iterate over dimensions
        
        int start = gridpoint[dimensionIndex] - 1;
        
        int end = gridpoint[dimensionIndex] + 1;
        
        for ( int index = start; index <= end; index++ )
        {
            gridpoint[dimensionIndex] = index;
            
            // guard against index error
            
            if (index < 0 || index > (dim.getSize()-1))
            {
                continue;
            }

            if ( dimensionIndex < (gridpoint.size()-1) )
            {
                // recurse
                find_neighbours(fs,arrayIndex,list,dimensionIndex+1,gridpoint);
            }
            else
            {
                // collect
                
                //std::cout << gridpoint << endl;
                
                typename Point<T>::ptr p = arrayIndex.get(gridpoint);
                
                if (p != NULL)
                {
                    list.push_back(p);
                }
            }
        }
        
        gridpoint[dimensionIndex] = start+1;
    }

    template <typename T>
    typename Point<T>::list
    ClusterList<T>::find_neighbours(FeatureSpace<T> *fs,
                                    ArrayIndex<T> &index,
                                    const typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        typename Point<T>::list neighbours;
        
        typename CoordinateSystem<T>::GridPoint gp = gridpoint;
        
        //cout << "finding neighbours of " << gridpoint << " :" << endl;
        
        this->find_neighbours(fs,index,neighbours,0,gp);
        
        return neighbours;
    }

    template <typename T>
    void
    ClusterList<T>::aggregate_zeroshifts(FeatureSpace<T> *fs,
                                         ArrayIndex<T> &index,
                                         bool show_progress)
    {
        boost::progress_display *progress = NULL;
        
        if (show_progress)
        {
            cout << endl << "Clustering zero-shift areas ...";
            progress = new boost::progress_display( fs->points.size() );
            start_timer();
        }

        for ( size_t i = 0; i < fs->points.size(); i++ )
        {
            if (show_progress)
            {
                progress->operator++();
            }
            
            M3DPoint<T> *current_point = (M3DPoint<T> *) fs->points[i];
            
            if (current_point->cluster != NULL)
            {
                continue;
            }
            
            // skip zeroshift to save time
            
            if (vector_norm(fs->spatial_component(current_point->shift)) == 0)
            {
                typename Point<T>::list zeroshift_neighbours;
                
                typename Point<T>::list neighbours = find_neighbours(fs,index,current_point->gridpoint);
                
                bool found_cluster = false;
                
                for (size_t ni = 0; ni < neighbours.size() && !found_cluster; ni++)
                {
                    M3DPoint<T> *n = (M3DPoint<T> *) neighbours.at(ni);
                    
                    if (n==current_point) continue;
                    
                    if (vector_norm(fs->spatial_component(n->shift)) == 0)
                    {
                        
                        if ( current_point->cluster == NULL && n->cluster == NULL )
                        {
                            // Neither current point nor neighbour have cluster
                            // => create new cluster
                            
                            size_t spatial_dims = fs->coordinate_system->size();
                            
                            typename Cluster<T>::ptr c = new Cluster<T>(current_point->values,spatial_dims);
                            
                            c->add_point(current_point);
                            
                            c->add_point(n);
                            
                            clusters.push_back(c);
                        }
                        else if ( current_point->cluster == NULL && n->cluster != NULL )
                        {
                            // neighbour has cluster
                            // => add current point to neighbour's cluster
                            n->cluster->add_point(current_point);
                        }
                        else if (current_point->cluster !=NULL && n->cluster == NULL)
                        {
                            // current point has cluster
                            // => add neighbour to current point's cluster
                            current_point->cluster->add_point(n);
                        }
                        else if ((current_point->cluster !=NULL && n->cluster != NULL)
                                 && (current_point->cluster != n->cluster ))
                        {
                            // current point's cluster and neighbour's cluster
                            // => merge current point's cluster into neighbour's cluster
                            typename Cluster<T>::ptr c = current_point->cluster;
                            
                            n->cluster->add_points(c->points,false);
                            
                            clusters.erase(find(clusters.begin(),clusters.end(),c));
                            
                            delete c;
                        }
                    }
                }
            }
        }
        
        for (size_t i=0; i < clusters.size(); i++)
        {
            typename Cluster<T>::ptr c = clusters.at(i);
            
            cout << "Found zeroshift cluster #" << i << " (" << c->points.size() << " points)." << endl;
        }
        
        if ( show_progress )
        {
            cout << "done. (Found " << clusters.size() << " clusters in " << stop_timer() << "s)" << endl;
            delete progress;
        }
    }
    
    template <typename T>
    void
    ClusterList<T>::check_clusters(FeatureSpace<T> *fs, ArrayIndex<T> &index)
    {
        // sanity checking
        
        typename Cluster<T>::list::iterator ci;
        
        for (ci=clusters.begin(); ci!=clusters.end(); ci++)
        {
            typename Cluster<T>::ptr c1 = *ci;
            
            typename Point<T>::list::iterator pi;
            
            for (pi = c1->points.begin(); pi != c1->points.end(); pi++)
            {
                M3DPoint<T> *p = (M3DPoint<T> *) *pi;
                
                M3DPoint<T> *pred = (M3DPoint<T> *) predecessor_of(fs,index,p);
                
                if (pred==NULL) continue;
                
                if (pred->cluster != c1)
                {
                    cout << "Point " << p->coordinate << " points to " << pred->coordinate << endl;
                    cout << "Their clusters are not identical!";
                    
//                    pred->cluster->add_points(c1->points,false);
//                    
//                    clusters.erase(ci);
//                    
//                    delete c1;
//                    
//                    // start over
//                    
//                    ci = clusters.begin();
                
                    typename Point<T>::list::iterator pci;
                    
                    for (pci = pred->cluster->points.begin(); pci != pred->cluster->points.end(); pci++)
                    {
                        M3DPoint<T> *pc = (M3DPoint<T> *) *pci;
                        
                        if (p->coordinate == pc->coordinate)
                        {
                            cout << "Point " << p << " @ " << p->coordinate
                                 << " has same coordinate as " << pc->coordinate << " @" << pc->coordinate
                                 << " but their clusters are not identical"
                                 << endl;
                        }
                    }
                }
            }
        }
    }
    
    template <typename T>
    void
    ClusterList<T>::aggregate_cluster_graph(const WeightFunction<T> *weight_function, FeatureSpace<T> *fs, const vector<T> &resolution, const bool& show_progress )
    {
        // PointIndex<T>::write_index_searches = true;
        
        boost::progress_display *progress = NULL;
        
        ArrayIndex<T> index(fs->coordinate_system,fs->points,false);
        
        this->aggregate_zeroshifts(fs,index,show_progress);
        
        // Sanity checking
        // this->check_clusters(fs,index);
        
        size_t cluster_id = this->clusters.size();
        
        if (show_progress)
        {
            cout << endl << "Analysing meanshift vector graph ...";
            start_timer();
            progress = new boost::progress_display( fs->points.size() );
        }
        
        for ( size_t i = 0; i < fs->points.size(); i++ )
        {
            if (show_progress)
            {
                progress->operator++();
            }
            
            M3DPoint<T> *current_point = (M3DPoint<T> *) fs->points[i];
            
            // skip zeroshift to save time
            
            if (vector_norm(fs->spatial_component(current_point->shift)) == 0)
            {
                continue;
            }
            
            // Find the predecessor through gridded shift
            
            size_t spatial_dims = fs->coordinate_system->size();
            
            vector<size_t> gridpoint(spatial_dims,0);
            
            for (size_t k=0; k < spatial_dims; k++)
            {
                gridpoint[k] = current_point->gridpoint[k] + current_point->gridded_shift[k];
            }
            
            M3DPoint<T> *predecessor = (M3DPoint<T> *)index.get(gridpoint);

            // Start testing
            
            if (predecessor != NULL)
            {
                #if DEBUG_GRAPH_AGGREGATION
                    cout << endl;
                    cout << "current point : " << current_point << " @ " << current_point->gridpoint << " (" << current_point->cluster << ")" << endl;
                    cout << "predecessor   : " << predecessor << " @ " << predecessor->gridpoint << " (" << predecessor->cluster << ")" << endl;
                    // cout << "(reverse lookup of " << x << " = " << gp << ")" << endl;
                #endif
                
                if (current_point->cluster == NULL && predecessor->cluster == NULL)
                {
                    // Neither point has a cluster
                    // => create new one
                    
                    typename Cluster<T>::ptr c = new Cluster<T>(current_point->values,fs->coordinate_system->size());
                    c->id = cluster_id++;
                    c->add_point(current_point);
                    c->add_point(predecessor);
                    clusters.push_back(c);
                    
                    #if DEBUG_GRAPH_AGGREGATION
                        cout << "created new cluster " << c << " (" << c->points.size() << " points)" << endl;
                    #endif
                }
                else if (current_point->cluster == NULL && predecessor->cluster != NULL)
                {
                    // current point has no cluster, but predecessor has one
                    // => add current point to predecessor's cluster
                    
                    predecessor->cluster->add_point(current_point);
                    #if DEBUG_GRAPH_AGGREGATION
                        cout << "added current point to cluster " << predecessor->cluster
                             << " (" << predecessor->cluster->points.size() << " points)" << endl;
                    #endif
                }
                else if (current_point->cluster != NULL && predecessor->cluster == NULL)
                {
                    // current point has a cluster, but predecessor has none
                    // => add predecessor to current point's cluster
                    
                    current_point->cluster->add_point(predecessor);

                    #if DEBUG_GRAPH_AGGREGATION
                    cout << "added predecessor to cluster " << current_point->cluster << " (" << current_point->cluster->points.size() << " points)" << endl;
                    #endif
                }
                else if (current_point->cluster != NULL
                         && predecessor->cluster != NULL
                         && (current_point->cluster != predecessor->cluster))
                {
                    // both points have different clusters
                    // => merge current cluster's points to predecessor's cluster
                    //    and delete current cluster
                    
                    typename Cluster<T>::ptr c1 = current_point->cluster;
                    
                    #if DEBUG_GRAPH_AGGREGATION
                        cout << "merging clusters " << c1 << " (" << c1->points.size() << " points)"
                             << "into " << predecessor->cluster << " (" << predecessor->cluster->points.size() << " points)"
                             << endl;
                    #endif
                    
                    // absorb predecessor
                
                    predecessor->cluster->add_points(c1->points,false);
                    
                    // remove it
                    
                    typename Cluster<T>::list::iterator fi;
                    
                    clusters.erase(find(clusters.begin(),clusters.end(),c1));
                    
                    delete c1;
                }
                else if (current_point->cluster != NULL
                         && predecessor->cluster != NULL
                         && (current_point->cluster == predecessor->cluster))
                {
                    // both points are already part of the same cluster
                    // => do nothing
                    
                    #if DEBUG_GRAPH_AGGREGATION
                        cout << "Both points are part of the same cluster. Skip." << endl;
                    #endif
                }
            }
        }
        
        // Iterate over the clusters and make them modes
        // the arithmetic mean
        
        for (size_t i=0; i < clusters.size(); i++)
        {
            typename Cluster<T>::ptr c = clusters.at(i);
            
            vector<T> mode = vector<T>( fs->feature_variables().size(), 0.0);
            
            for (size_t j=0; j < c->points.size(); j++)
            {
                M3DPoint<T> *p = (M3DPoint<T> *) c->points.at(j);
                
                mode += p->values;
            }
            
            mode /= ((T) c->points.size());
            
            c->mode = mode;
            
            c->id = i;
            
            cout << "Found zeroshift cluster #" << i << " at " << mode << " (" << c->points.size() << " points)." << endl;
        }
        
        // Sanity checking
        // this->check_clusters(fs,index);
        
        if ( show_progress )
        {
            cout << "done. (Found " << clusters.size() << " clusters in " << stop_timer() << "s)" << endl;
            delete progress;
        }
        
        // PointIndex<T>::write_index_searches = false;
    }
    
    template <typename T>
    typename Cluster<T>::list
    ClusterList<T>::neighbours_of( typename Cluster<T>::ptr cluster, PointIndex<T> *index, const vector<T> &resolution )
    {
        typename Cluster<T>::list neighbouring_clusters;
        
        RangeSearchParams<T> search_params( resolution );
        
        typename Point<T>::list::const_iterator pi;
        
        for ( pi = cluster->points.begin(); pi != cluster->points.end(); pi++ )
        {
        	typename M3DPoint<T>::ptr p = (M3DPoint<T> *) *pi;
            
            typename Point<T>::list *neighbours = index->search( p->values, &search_params );
            
            typename Point<T>::list::const_iterator ni;
            
            for ( ni = neighbours->begin(); ni != neighbours->end(); ni++ )
            {
                M3DPoint<T> *n = (M3DPoint<T> *) *ni;
                
                // Exclude points, that have not been clustered.
                // This can happen because scale-space filtering
                // creates new points, but those are not associated
                // with clusters in later steps
                if ( n->cluster == NULL ) continue;
                
                if ( n->cluster != p->cluster )
                {
                    typename Cluster<T>::list::const_iterator fi = find( neighbouring_clusters.begin(), neighbouring_clusters.end(), n->cluster );
                    
                    if ( fi == neighbouring_clusters.end() )
                    {
                        neighbouring_clusters.push_back( n->cluster );
                    }
                }
            }
            
            delete neighbours;
        }
        
        return neighbouring_clusters;
    }
    
    template <typename T>
    void
    ClusterList<T>::get_boundary_points(typename Cluster<T>::ptr c1,
                                        typename Cluster<T>::ptr c2,
                                        typename Point<T>::list &boundary_points,
                                        PointIndex<T> *index,
                                        const vector<T> &resolution )
    {
        RangeSearchParams<T> search_params( m_neighbourhood_range_search_multiplier * resolution );
        
        typename Point<T>::list::const_iterator pi;
        
        for ( pi = c1->points.begin(); pi != c1->points.end(); pi++ )
        {
            typename Point<T>::ptr p = *pi;
            
            typename Point<T>::list *neighbours = index->search( p->values, &search_params );
            
            typename Point<T>::list::const_iterator ni;
            
            for ( ni = neighbours->begin(); ni != neighbours->end(); ni++ )
            {
                typename M3DPoint<T>::ptr n = (M3DPoint<T> *) *ni;
                
                if ( n->cluster == c2 )
                {
                    // check every time to avoid double adding
                    
                    typename Point<T>::list::const_iterator fi = find( boundary_points.begin(), boundary_points.end(), n );
                    
                    if ( fi == boundary_points.end() )
                    {
                        boundary_points.push_back( n );
                    }
                    
                    fi = find( boundary_points.begin(), boundary_points.end(), p );
                    
                    if ( fi == boundary_points.end() )
                    {
                        boundary_points.push_back( p );
                    }
                }
            }
            
            delete neighbours;
        }

        for ( pi = c2->points.begin(); pi != c2->points.end(); pi++ )
        {
            typename Point<T>::ptr p = *pi;
            
            typename Point<T>::list *neighbours = index->search( p->values, &search_params );
            
            typename Point<T>::list::const_iterator ni;
            
            for ( ni = neighbours->begin(); ni != neighbours->end(); ni++ )
            {
            	typename M3DPoint<T>::ptr n = (M3DPoint<T> *) *ni;
                
                if ( n->cluster == c1 )
                {
                    // check every time to avoid double adding
                    
                    typename Point<T>::list::const_iterator fi = find( boundary_points.begin(), boundary_points.end(), n );
                    
                    if ( fi == boundary_points.end() )
                    {
                        boundary_points.push_back( n );
                    }
                    
                    fi = find( boundary_points.begin(), boundary_points.end(), p );
                    
                    if ( fi == boundary_points.end() )
                    {
                        boundary_points.push_back( p );
                    }
                }
            }
            
            delete neighbours;
        }
    };
    
    // Couple of macros
    
    template <typename T>
    bool
    ClusterList<T>::should_merge_neighbouring_clusters( typename Cluster<T>::ptr c1,
                                                        typename Cluster<T>::ptr c2,
                                                        const WeightFunction<T> *weight_function,
                                                        PointIndex<T> *index,
                                                        const vector<T> &resolution,
                                                        const double &drf_threshold )
    {
        bool should_merge = false;
        
        // Collect the points on each side of the boundary
        
        typename Point<T>::list boundary_points;
        
        this->get_boundary_points( c1, c2, boundary_points, index, resolution );
        
        // if no common boundary exists, the point is moot
        
#if DEBUG_CLUSTER_MERGING_DECISION
        std::cout << "===================================================================" << std::endl;
        std::cout << "comparing #" << c1->id << " = " << c1->mode << " (" << c1->size() << " points)" << std::endl;
        std::cout << "with      #" << c2->id << " = " << c2->mode << " (" << c2->size() << " points)" << std::endl;
        
#endif
        if ( boundary_points.size() > 0 )
        {
            // the classification of dynamic range of the signal on each side
            // of the boundary
            
            T c1_dyn_range_factor = dynamic_range_factor( c1, boundary_points, weight_function );
            
            T c2_dyn_range_factor = dynamic_range_factor( c2, boundary_points, weight_function );
            
            bool dynamic_range_test = sup( c1_dyn_range_factor, c2_dyn_range_factor ) >= drf_threshold;
            
#if DEBUG_CLUSTER_MERGING_DECISION
            cout << "\tc1-drf = " << c1_dyn_range_factor << " , c2-drf = " << c2_dyn_range_factor << endl;
#endif
            
            // Indicators for merging are: high dynamic range on both sides of the boundary
            // as well as high signal variablility in both

            should_merge =  dynamic_range_test;
        }
#if DEBUG_CLUSTER_MERGING_DECISION
        else
        {
            cout << "no common boundary" << endl;
        }
#endif 
        
#if DEBUG_CLUSTER_MERGING_DECISION
        std::cout << "\t==> should merge = " << (should_merge ? "yes" : "no") << std::endl;
        std::cout << "===================================================================" << std::endl;
#endif

        return should_merge;
    }
    
    template <typename T>
    void
    ClusterList<T>::write_boundaries(const WeightFunction<T> *weight_function,
                                     FeatureSpace<T> *fs,
                                     PointIndex<T> *index,
                                     const vector<T> &resolution )
    {
        // collate the data
        
        typedef vector< typename Point<T>::list > boundaries_t;
        
        boundaries_t boundaries;
        
        
        typedef vector< std::string > boundary_key_t;
        
        boundary_key_t boundary_keys;
        
        
        vector<T> var_c1, var_c2, var_boundary;
        
        vector<T> range_factor_c1, range_factor_c2;
        
        vector< typename Cluster<T>::id_t >  cluster_index_1, cluster_index_2;
        
        typename Cluster<T>::list::const_iterator ci;

        for ( ci = clusters.begin(); ci != clusters.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            typename Cluster<T>::list neighbours = neighbours_of( c, index, resolution, weight_function );
            
            if ( neighbours.size() > 0 )
            {
                // go over the list of neighbours and find candidates for merging
                
                typename Cluster<T>::list::const_iterator ni;
                
                for ( ni = neighbours.begin(); ni != neighbours.end(); ni++ )
                {
                    typename Cluster<T>::ptr n = *ni;

                    std::string key = boost::lexical_cast<string>(inf(c->id,n->id)) + "-" + boost::lexical_cast<string>(sup(c->id,n->id));
                    
                    typename boundary_key_t::const_iterator fi = find( boundary_keys.begin(), boundary_keys.end(), key );
                    
                    if ( fi == boundary_keys.end() )
                    {
                        boundary_keys.push_back( key );
                        
                        typename Point<T>::list boundary_points;
                        
                        this->get_boundary_points( c, n, boundary_points, index, resolution );
                        
                        if ( boundary_points.size() == 0 ) continue;
                            
                        
                        boundaries.push_back( boundary_points );
                        

                        var_boundary.push_back( relative_variability( weight_function, boundary_points ) );
                        
                        var_c1.push_back( relative_variability( weight_function, c->points ) );

                        var_c2.push_back( relative_variability( weight_function, n->points ) );

                        
                        range_factor_c1.push_back( dynamic_range_factor( c, boundary_points, weight_function ) );

                        range_factor_c2.push_back( dynamic_range_factor( n, boundary_points, weight_function ) );
                        
                        
                        cluster_index_1.push_back(c->id);
                        
                        cluster_index_2.push_back(n->id);
                    }
                }
            }
        }
        
        for ( size_t index = 0; index < boundaries.size(); index++ )
        {
            typename Point<T>::list b = boundaries[index];
            std::string fn = fs->filename() + "_boundary_" + boost::lexical_cast<string>(index) + ".vtk";
            boost::replace_all( fn, "/", "_" );
            boost::replace_all( fn, "..", "" );
            cfa::utils::VisitUtils<T>::write_pointlist_vtk(fn, &b, fs->coordinate_system->size());
        }
    
        std::string fn = fs->filename() + "_boundary_correlations.txt";
        std::ofstream f( fn.c_str() );
        f << "#\t"
        << "c1\t"
        << "c2\t"
//        << "var_c1\t"
//        << "var_c2\t"
//        << "var_b\t"
        << "drf_1\t"
        << "drf_2\t"
	  << std::endl;
        for ( size_t index = 0; index < boundaries.size(); index++ )
        {
            f << index << "\t"
            << cluster_index_1[index] << "\t"
            << cluster_index_2[index] << "\t"
//            << var_c1[index] << "\t"
//            << var_c2[index] << "\t"
//            << var_boundary[index] << "\t"
            << range_factor_c1[index] << "\t"
	      << range_factor_c2[index] << std::endl;
        }
    }

    
    template <typename T>
    T
    ClusterList<T>::relative_variability( const WeightFunction<T> *weight_function, const typename Point<T>::list &points )
    {
        // calculate mean
        
        T mean = 0.0;
        
        typename Point<T>::list::const_iterator pi;
        
        for ( pi = points.begin(); pi != points.end(); pi++ )
        {
            typename Point<T>::ptr p = *pi;
            
            T value = weight_function->operator()(p->coordinate,p->values);
            
            mean += value;
        }
        
        mean /= ((T)points.size());
        
        // calculate standard deviation
        
        T standard_deviation = 0.0;
        
        for ( pi = points.begin(); pi != points.end(); pi++ )
        {
            typename Point<T>::ptr p = *pi;
            
            T value = weight_function->operator()(p->coordinate,p->values);
        
            standard_deviation += ( value - mean ) * ( value - mean );
        }
        
        standard_deviation = sqrt( standard_deviation / ((T)points.size()) );
        
        // relative deviation
        
        T cv = standard_deviation / mean;
        
#if DEBUG_CLUSTER_MERGING
        std::cout << "Sample of " << points.size() << " points in variable index " << variable_index
             << " has m=" << mean << ", s=" << standard_deviation << ", cv=" << cv << std::endl;
#endif
        return cv;
    }

    template <typename T>
    T
    ClusterList<T>::dynamic_range_factor( typename Cluster<T>::ptr cluster,
    									  const typename Point<T>::list &points,
    									  const WeightFunction<T> *weight_function )
    {
        // Obtain the dynamic range if the list
        
        T points_lower_bound, points_upper_bound;
        
        Cluster<T>::dynamic_range( points, weight_function, points_lower_bound, points_upper_bound );
        
        T cluster_lower_bound, cluster_upper_bound;
        
        cluster->dynamic_range( weight_function, cluster_lower_bound, cluster_upper_bound );

        // Compare the two ranges
        
//        T p_range = points_upper_bound - points_lower_bound;
//        
//        T c_range = cluster_upper_bound - cluster_lower_bound;
//        
//        T factor = ( p_range == 0 || c_range == 0 ) ? 0.0 : inf(p_range,c_range) / sup(p_range,c_range);

        // Correction: compare upper values only and see if they're in the same
        // ballpark
        
        T factor = inf( points_upper_bound, cluster_upper_bound ) / sup( points_upper_bound, cluster_upper_bound );
        
#if DEBUG_CLUSTER_MERGING
        std::cout << "-- Dynamic Range Comparison " << std::endl;
        std::cout << "\tcluster at " << cluster->mode << "(" << cluster->size() << " points) compared to boundary line (" << points.size() << " points)" << std::endl;
        std::cout << "\t\tcluster  dynamic range = [" << cluster_lower_bound << "," << cluster_upper_bound << "]" << std::endl;
        std::cout << "\t\tboundary dynamic range = [" << points_lower_bound << "," << points_upper_bound << "]" << std::endl;
        std::cout << "\tfactor = " << factor << std::endl;
#endif
        return factor;
    }

    template <typename T>
    typename Cluster<T>::ptr
    ClusterList<T>::merge_clusters( typename Cluster<T>::ptr c1, typename Cluster<T>::ptr c2 )
    {
        vector<T> merged_mode = (T)0.5 * ( c1->mode + c2->mode );
        
        typename Cluster<T>::ptr merged_cluster = new Cluster<T>( merged_mode, this->dimensions.size() );
        
        merged_cluster->add_points( c1->points );
        
        merged_cluster->add_points( c2->points );
        
        merged_cluster->id = c1->id;
        
        if (c1->m_weight_range_calculated && c2->m_weight_range_calculated)
        {
            merged_cluster->m_min_weight = inf(c1->m_min_weight, c2->m_min_weight);
            merged_cluster->m_max_weight = sup(c1->m_max_weight, c2->m_max_weight);
            merged_cluster->m_weight_range_calculated = true;
        }
        
        return merged_cluster;
    }

    
    template <typename T>
    void
    ClusterList<T>::aggregate_clusters_by_boundary_analysis( const WeightFunction<T> *weight_function,
                                                             PointIndex<T> *index,
                                                             const vector<T> &resolution,
                                                             const double &drf_threshold,
                                                             const bool &show_progress )
    {
        // No merging
        if (drf_threshold==1.0) return;
        
        // Provide a second index for the merge examination. They use a smaller
        // resolution and would trigger re-calculation of the index every other
        // step. It's cheaper to provide them with their own index
        
        PointIndex<T> *index2 = PointIndex<T>::create( (FeatureSpace<T> *) index->feature_space() );
        
        ConsoleSpinner *spinner = NULL;

        if ( show_progress )
        {
            cout << endl << "Aggregating clusters into objects ..." << endl;
            spinner = new ConsoleSpinner();
            start_timer();
        }
        
        bool had_merge = false;
        
        size_t merge_count = 0;
        

        typedef vector< std::string > boundary_key_t;
        
        boundary_key_t boundary_keys;
        
        typename Cluster<T>::list scheduledForRemoval;
        
        size_t pass_counter = 1;
        
        do
        {
            if ( show_progress )
            {
                cout << "\tpass #" << pass_counter << " : " << clusters.size() << " -> ";
            }

            had_merge = false;
            
            // Start at the top and work your way down
            
            typename Cluster<T>::list::const_iterator ci;
            
            for ( ci = clusters.begin(); ci != clusters.end() && ! had_merge; ci++ )
            {
                typename Cluster<T>::ptr c = *ci;
                
                typename Cluster<T>::list neighbours = neighbours_of( c, index, resolution );
                
                Cluster<T> *merged = c;

                // go over the list of neighbours and find candidates for merging
                // Merge all suitable candidates of c's neighbours together
                
                typename Cluster<T>::list::iterator ni;
                
                for ( ni = neighbours.begin(); ni != neighbours.end(); ni++ )
                {
                    typename Cluster<T>::ptr n = *ni;
                    
                    if (should_merge_neighbouring_clusters( c, n, weight_function, index2, resolution, drf_threshold))
                    {
                        scheduledForRemoval.push_back(n);
                        
                        merged = merge_clusters( merged, n );
                    
                        had_merge = true;
                    
                        merge_count++;
                    }
                }
                
                if (had_merge)
                {
                    scheduledForRemoval.push_back(c);
                }
                
                // remove c and all merged neighbours from the list,
                // replace them by the merged one and do it over again
                // until no more merges happened
                
                if (had_merge)
                {
                    for (size_t i=0; i < scheduledForRemoval.size(); i++)
                    {
                        Cluster<T> *rc = scheduledForRemoval[i];
                        typename Cluster<T>::list::iterator fi = find( clusters.begin(), clusters.end(), rc );
                        if (fi!=clusters.end())
                        {
                            clusters.erase( fi );
                        }
                    }
                    
                    clusters.push_back(merged);
                    
                    break;
                }
            }
            
            if ( show_progress )
            {
                cout << clusters.size() << endl;
            }
            
            pass_counter++;
            
        } while (had_merge);
        
        
        // Delete the removed clusters
        
        // remove duplicates
        
        sort( scheduledForRemoval.begin(), scheduledForRemoval.end() );
        scheduledForRemoval.erase( unique( scheduledForRemoval.begin(), scheduledForRemoval.end() ), scheduledForRemoval.end() );
        
        //    delete

        for (size_t i=0; i < scheduledForRemoval.size(); i++)
        {
            Cluster<T> *rc = scheduledForRemoval[i];
            delete rc;
        }
        
        if ( show_progress )
        {
            delete spinner;
            cout << "done. (Found " << clusters.size() << " objects in " << stop_timer() << "s, " << merge_count << " merges)" << endl;
        }
        
        delete index2;
    };
    
    template <typename T>
    void
    ClusterList<T>::erase_identifiers()
    {
        for ( size_t i=0; i < clusters.size(); i++ )
        {
            clusters[i]->id = Cluster<T>::NO_ID;
        }
    };
    
    template <typename T>
    void
    ClusterList<T>::retag_identifiers()
    {
        for ( size_t i=0; i < clusters.size(); i++ )
        {
            clusters[i]->id = i;
        }
    };
    
    
    template <typename T>
    void
    ClusterList<T>::reset_clustering( FeatureSpace<T> *fs )
    {
        struct clear_cluster
        {
            void operator() (void *p)
            {
                static_cast< M3DPoint<T> * >(p)->cluster = NULL;
            };
        } clear_cluster;
        
        for_each( fs->points.begin(), fs->points.end(), clear_cluster );
    };
    
    template <typename T>
    void
    ClusterList<T>::sanity_check( const FeatureSpace<T> *fs )
    {
        size_t point_count = 0;
        
        for ( size_t i=0; i < clusters.size(); i++ )
        {
            point_count += clusters[i]->points.size();
        }
        
        assert( point_count == fs->size() );
    };
    
#pragma mark -
#pragma mark Coalescence Merging
    
    template <typename T>
    bool
    ClusterList<T>::are_neighbours(const Cluster<T> *c1,
                                   const Cluster<T> *c2,
                                   PointIndex<T> *index,
                                   const vector<T> &resolution)
    {
        bool isNeighbour = false;
        
        RangeSearchParams<T> search_params( resolution );
        
        typename Point<T>::list::const_iterator pi;
        
        for ( pi = c1->points.begin(); pi != c1->points.end(); pi++ )
        {
        	typename M3DPoint<T>::ptr p = (M3DPoint<T> *) *pi;
            
            typename Point<T>::list *neighbours = index->search( p->values, &search_params );
            
            typename Point<T>::list::const_iterator ni;
            
            for ( ni = neighbours->begin(); ni != neighbours->end(); ni++ )
            {
                M3DPoint<T> *n = (M3DPoint<T> *) *ni;
                
                if ( n->cluster == c2 )
                {
                    isNeighbour = true;
                    break;
                }
            }
        }
        
        return isNeighbour;
    }

    
    // Sort all clusters in ascending order by weight response

    template <typename T>
    class ModalWeightComparator
    {
    private:
        const WeightFunction<T> *m_weight;
    public:
        ModalWeightComparator(const WeightFunction<T> *w) {m_weight=w;}
        
        bool
        operator() (const Cluster<T> *c1, const Cluster<T> *c2)
        {
            T w1 = c1->modal_weight_response(m_weight);
            T w2 = c2->modal_weight_response(m_weight);
            return w1 < w2;
        }
    };
    
    
    template <class T>
    void
    ClusterList<T>::aggregate_by_coalescence(const WeightFunction<T> *weight_function,
                                             PointIndex<T> *index,
                                             const vector<T> &resolution,
                                             size_t passes,
                                             bool show_progress)
    {
        
        // Provide a second index for the merge examination. They use a smaller
        // resolution and would trigger re-calculation of the index every other
        // step. It's cheaper to provide them with their own index
        
        PointIndex<T> *index2 = PointIndex<T>::create( (FeatureSpace<T> *) index->feature_space() );

        // Sort by modal weight function ascending
        
        ModalWeightComparator<T> comp(weight_function);
        
        sort(clusters.begin(), clusters.end(), comp);
        
        // Traverse the list and find the closest neighbour
        // with higher weight response.
        
        typedef Cluster<T> *            cptr_t;
        typedef set<cptr_t>             cptr_set_t;
        typedef map<cptr_t,cptr_t>      cptr_match_t;
        typedef map<cptr_t,cptr_set_t>  cptr_map_t;
        
        
        // Coalescence
        
        for (int pass=0; pass < passes; pass++)
        {
            // pre-calculate direct neighbours
            
            boost::progress_display *progress = NULL;
            
            if (show_progress)
            {
                cout << endl << "Analysing neighbourhood relationships ...";
                start_timer();
                progress = new boost::progress_display( clusters.size() );
            }
            
            cptr_map_t neighbours;
            
            for (size_t ci=0; ci < (clusters.size()-1); ci++)
            {
                Cluster<T> *c1 = clusters[ci];
                
                if (show_progress)
                {
                    progress->operator++();
                }
                
                typename Cluster<T>::list direct_neighbours;
                direct_neighbours = neighbours_of( c1, index2, resolution );
                
                if ( direct_neighbours.size() > 0 )
                {
                    // go over the list of neighbours and find candidates for merging
                    
                    typename Cluster<T>::list::const_iterator ni;
                    
                    for ( ni = direct_neighbours.begin(); ni != direct_neighbours.end(); ni++ )
                    {
                        cptr_t n = *ni;

                        neighbours[c1].insert(n);
                        
                        neighbours[n].insert(c1);
                    }
                }
            }
            
            if (show_progress)
            {
                cout << endl << "done. " << stop_timer() << "s)" << endl;
                delete progress;
            }

            // perform the coalescence search
            
            cptr_match_t matches;

            if (show_progress)
            {
                cout << endl << "Finding matches ...";
                start_timer();
                progress = new boost::progress_display( clusters.size());
            }
            
            for (size_t ci=0; ci < (clusters.size()-1); ci++)
            {
                if (show_progress)
                {
                    progress->operator++();
                }

                Cluster<T> *c1 = clusters[ci];
                
                // search 'upwards'
                
                int highest_index = -1;
                
                for (int ni=ci+1; ni < clusters.size(); ni++)
                {
                    cptr_t c2 = clusters[ni];
                    
                    typename cptr_set_t::iterator fi;
                    fi = find( neighbours[c1].begin(), neighbours[c1].end(), c2);
                    
                    if ( fi != neighbours[c1].end() )
                    {
                        // found an "upstream" neighbour
                        highest_index = ni;
                    }
                }
                
                // And the winner is ...
                
                if (highest_index > 0)
                {
                    matches[c1] = clusters[highest_index];
                }
            }
            
            if (show_progress)
            {
                cout << endl << "done. " << stop_timer() << "s)" << endl;
                delete progress;
            }
        
            // coalesce
            
            if (show_progress)
            {
                cout << endl << "Coalesce clusters ...";
                start_timer();
                progress = new boost::progress_display( matches.size() );
            }
            
            cptr_match_t merged_map;
            
            vector<cptr_t> used;
            
            typename cptr_match_t::iterator mi;
            
#if DEBUG_COALESCENCE_MATCHING
            for (mi = matches.begin(); mi != matches.end(); mi++)
            {
                cout << "Cluster #" << mi->first->id << " goes to " << mi->second->id << endl;
            }
#endif
            
            for (mi = matches.begin(); mi != matches.end(); mi++)
            {
                if (show_progress)
                {
                    progress->operator++();
                }

                //
                // Figure out, which one to use for merging
                //
                
                // Was this cluster merged with another cluster
                // before?
                
                // TODO: this should not be necessary, because
                // the clusters are ordered and the upper cluster
                // never comes before the lower
                
                typename cptr_match_t::iterator lc_find;
                lc_find = merged_map.find(mi->first);

                // if lower cluster was previously merged, use
                // the merged one. Else use the lower cluster
                // itself
                
                cptr_t lower_cluster = (lc_find == merged_map.end())
                                     ? mi->first
                                     : lc_find->second;
                
                //
                // Now for the 'matched' one.
                //
                
                cptr_t upper_cluster = mi->second;
                    
                // Was the upper cluster merged with another
                // cluster before?
                
                typename cptr_match_t::const_iterator uc_find;
                uc_find = merged_map.find(upper_cluster);

                if (uc_find != merged_map.end())
                {
                    // Hold on to the previous merged
                    cptr_t previous_m = uc_find->second;
                    
                    // merge with the result of the
                    // previous merge
                    cptr_t m = merge_clusters(lower_cluster,previous_m);
                    
                    // replace all occurences of previous merge
                    // in the merged map with the new one
                    
                    typename cptr_match_t::iterator it;
                    
                    for (it=merged_map.begin(); it!=merged_map.end(); it++)
                    {
                        if (it->second == previous_m)
                        {
                            it->second = m;
                        }
                    }
                    
                    // dispose of the previous matched
                    
                    delete previous_m;
                }
                else
                {
                    // merge
                    cptr_t m = merge_clusters(lower_cluster,upper_cluster);
                    
                    // create a new entries in the matching map
                    merged_map[lower_cluster] = m;
                    merged_map[upper_cluster] = m;
                }
            }
            
            if (show_progress)
            {
                cout << endl << "done. " << stop_timer() << "s)" << endl;
                delete progress;
            }
            
            // Extract the merged results by uniquing the merged
            // map's value set
            
            vector<cptr_t> result;
            
            for (typename cptr_match_t::iterator it=merged_map.begin(); it!=merged_map.end(); it++)
            {
                result.push_back(it->second);
            }
            
            // remove duplicates
            
            sort(result.begin(),result.end());
            unique(result.begin(),result.end());
            
            // release original clusters and replace
            // with the result
            
            while (!clusters.empty())
            {
                cptr_t c = clusters.back();
                clusters.pop_back();
                delete c;
            }
            
            clusters = result;
            
        }
    }

    


}; //namespace

#endif
