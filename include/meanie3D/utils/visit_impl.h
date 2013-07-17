#ifndef _M3D_VisualizationUtils_Impl_H_
#define _M3D_VisualizationUtils_Impl_H_

#include <iostream>
#include <stdlib.h>
#include <vector>

#include <meanie3D/types/cluster.h>
#include <meanie3D/types/cluster_list.h>

namespace m3D { namespace utils {

	using namespace std;

    template <typename T>
	vector<size_t>
    VisitUtils<T>::VTK_DIMENSION_INDEXES;

	/** Write out the cluster list as list of points. The point's 'value' is the size of the cluster
     * @param full path to filename, including extension '.vtk'
     * @param cluster list
     */
    template <typename T>
    void
    VisitUtils<T>::write_cluster_modes_vtk( const string &filename, const typename Cluster<T>::list &list, bool spatial_only )
    {
        ofstream f( filename.c_str() );
        f << fixed << setprecision(4);

        // Write Header

        f << "# vtk DataFile Version 3.0" << endl;
        f << "point list " << endl;
        f << "ASCII" << endl;
        f << "DATASET UNSTRUCTURED_GRID" << endl;
        f << "POINTS " << list.size() << " FLOAT" << endl;
        
        for ( size_t ci = 0; ci < list.size(); ci++ )
        {
            typename Cluster<T>::ptr c = list[ci];
            
            size_t dim_count = spatial_only ? c->spatial_dimension() : c->dimension();
            
            vector<T> mode = c->mode;
            
            size_t dims_plotted = 0;

            for ( size_t vi = 0; vi < dim_count; vi++)
            {
                size_t dim_index = vi;
                
                if ( vi < VTK_DIMENSION_INDEXES.size() )
                {
                    dim_index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                }

                f << mode[dim_index] << "\t";
                
                dims_plotted++;
            }

            while (dims_plotted < 3)
            {
                f << "0.0\t";
                
                dims_plotted++;
            }

            f << endl;
        }

        f << endl;
        f << "POINT_DATA " << list.size() << endl;
        f << "SCALARS mode INT" << endl;
        f << "LOOKUP_TABLE default" << endl;
        for ( size_t pi = 0; pi < list.size() ; pi++ )
        {
            typename Cluster<T>::ptr c = list[pi];
            
            size_t dim_count = spatial_only ? c->spatial_dimension() : c->dimension();

            if ( spatial_only )
            {
                f << c->id << endl;
            }
            else
            {
                f << c->mode[ dim_count-1 ];
            }
        }

        f.close();
    }
    
    template <class T>
    void
    VisitUtils<T>::write_geometrical_cluster_centers_vtk(const string &filename,
                                                         const typename Cluster<T>::list &list)
    {
        ofstream f( filename.c_str() );
        f << fixed << setprecision(4);
        
        // determine the spatial dimension
        typename Cluster<T>::ptr c = list[0];
        typename Point<T>::ptr p = c->points[0];
        size_t spatial_dims = p->coordinate.size();

        // Write Header
        
        f << "# vtk DataFile Version 3.0" << endl;
        f << "point list " << endl;
        f << "ASCII" << endl;
        f << "DATASET UNSTRUCTURED_GRID" << endl;
        f << "POINTS " << list.size() << " FLOAT" << endl;
        
        for ( size_t ci = 0; ci < list.size(); ci++ )
        {
            typename Cluster<T>::ptr c = list[ci];
            
            // obtain the geometrical center
            
            vector<T> mode = c->geometrical_center(spatial_dims);
            
            size_t dims_plotted = 0;
            
            for ( size_t vi = 0; vi < spatial_dims; vi++)
            {
                size_t dim_index = vi;
                
                if ( vi < VTK_DIMENSION_INDEXES.size() )
                {
                    dim_index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                }
                
                f << mode[dim_index] << "\t";
                
                dims_plotted++;
            }
            
            while (dims_plotted < 3)
            {
                f << "0.0\t";
                
                dims_plotted++;
            }
            
            f << endl;
        }
        
        f << endl;
        f << "POINT_DATA " << list.size() << endl;
        f << "SCALARS geometrical_center INT" << endl;
        f << "LOOKUP_TABLE default" << endl;
        for ( size_t pi = 0; pi < list.size() ; pi++ )
        {
            typename Cluster<T>::ptr c = list[pi];
            
            f << c->id << endl;
        }
        
        f.close();
    }


    template <class T>
    void
    VisitUtils<T>::write_clusters_vtk(const string &base_name,
                                      const typename Cluster<T>::list &list,
                                      const vector<T> &bandwidths,
                                      bool use_ids,
                                      bool only_boundary)
    {
        string basename = base_name;

        for ( size_t ci = 0; ci < list.size(); ci++ )
        {
            vector<T> mode = list[ci]->mode;

            // Write cluster out

            boost::replace_all( basename, "/", "_" );
            boost::replace_all( basename, "..", "" );

            string filename = basename + "_cluster_" + boost::lexical_cast<string>( use_ids ? list[ci]->id : ci ) + ".vtk";

            ofstream f( filename.c_str() );
            f << fixed << setprecision(4);
            
            size_t num_points = list[ci]->points.size();
            
            // If it's boundary points only, we need to count them upfront
            
            if (only_boundary)
            {
                num_points = 0;
                
                for ( size_t pi = 0; pi < list[ci]->points.size(); pi++ )
                {
                    M3DPoint<T> *p = (M3DPoint<T> *) list[ci]->points[pi];
                    
                    if (p->isBoundary)
                    {
                        num_points++;
                    }
                }
            }

            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Meanshift Clustering Result" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << num_points << " FLOAT" << endl;

            // Write point coordinates out as unstructured grid

            size_t point_dim = list[ci]->points[0]->coordinate.size();

            for ( size_t pi = 0; pi < list[ci]->points.size(); pi++ )
            {
                M3DPoint<T> *p = (M3DPoint<T> *) list[ci]->points[pi];
                
                if (only_boundary && !p->isBoundary) continue;

                for ( size_t vi = 0; vi < point_dim; vi++)
                {
                    size_t index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];

                    f << p->coordinate[index] << "\t";
                }

                if ( point_dim < 3 )
                {
                    f << "0.0";
                }

                f << endl;
            }

            // Write point data out. Only take the first value after coordinates

            f << endl;
            f << "POINT_DATA " << num_points << endl;
            f << "SCALARS cluster FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;

            for ( size_t pi = 0; pi < list[ci]->points.size(); pi++ )
            {
                M3DPoint<T> *p = (M3DPoint<T> *) list[ci]->points[pi];
                
                if (only_boundary && !p->isBoundary) continue;

                // f << p->trajectory_length << endl;

                f << p->values[ point_dim ] << endl;
            }

            f.close();
        }
    };
    
    template <class T>
    void
    VisitUtils<T>::write_clusters_vtk( typename ClusterList<T>::ptr list, std::string infix, bool only_boundary )
    {
        string basename = list->source_file;
        
        for ( size_t ci = 0; ci < list->size(); ci++ )
        {
            vector<T> mode = list->clusters[ci]->mode;
            
            // Write cluster out
            
            boost::replace_all( basename, "/", "_" );
            boost::replace_all( basename, "..", "" );
            
            string filename = basename + infix + boost::lexical_cast<string>( list->clusters[ci]->id ) + ".vtk";
            
            ofstream f( filename.c_str() );
            f << fixed << setprecision(4);
            
            size_t num_points = list->clusters[ci]->points.size();
            
            // If it's boundary points only, we need to count them upfront
            
            if (only_boundary)
            {
                num_points = 0;
                
                for ( size_t pi = 0; pi < list->clusters[ci]->points.size(); pi++ )
                {
                    M3DPoint<T> *p = (M3DPoint<T> *) list->clusters[ci]->points[pi];
                    
                    if (p->isBoundary)
                    {
                        num_points++;
                    }
                }
            }
            
            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Meanshift Clustering Result" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << num_points << " FLOAT" << endl;
            
            // Write point coordinates out as unstructured grid
            
            size_t point_dim = list->clusters[ci]->points[0]->coordinate.size();
            
            for ( size_t pi = 0; pi < list->clusters[ci]->points.size(); pi++ )
            {
                M3DPoint<T> *p = (M3DPoint<T> *) list->clusters[ci]->points[pi];
                
                if (only_boundary && !p->isBoundary) continue;
                
                for ( size_t vi = 0; vi < point_dim; vi++)
                {
                    size_t index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                    
                    f << p->coordinate[index] << "\t";
                }
                
                if ( point_dim < 3 )
                {
                    f << "0.0";
                }
                
                f << endl;
            }
            
            // Write point data out. Only take the first value after coordinates
            
            f << endl;
            f << "POINT_DATA " << num_points << endl;
            f << "SCALARS cluster FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;
            
            for ( size_t pi = 0; pi < list->clusters[ci]->points.size() - 1 ; pi++ )
            {
                M3DPoint<T> *p = (M3DPoint<T> *) list->clusters[ci]->points[pi];
                
                if (only_boundary && !p->isBoundary) continue;
                
                // f << p->trajectory_length << endl;
                
                f << p->values[ point_dim ] << endl;
            }
            
            f.close();
        }
    };
    
    template <class T>
    void
    VisitUtils<T>::write_center_tracks_vtk(typename Tracking<T>::trackmap_t &track_map,
                                           const std::string &basename,
                                           size_t spatial_dimensions)
    {
        for (typename Tracking<T>::trackmap_t::iterator tmi = track_map.begin(); tmi != track_map.end(); tmi++)
        {
            typename Tracking<T>::track_t *track = tmi->second;
            
            string filename = basename + "-track_" + boost::lexical_cast<string>( tmi->first ) + ".vtk";
            
            ofstream f( filename.c_str() );
            f << fixed << setprecision(4);
            
            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Meanie3D Track" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << track->size() << " FLOAT" << endl;

            // Write point coordinates out as unstructured grid
            
            for ( typename Tracking<T>::track_t::iterator ti=track->begin(); ti!=track->end(); ti++ )
            {
                vector<T> center = ti->geometrical_center(spatial_dimensions);

                for ( size_t vi = 0; vi < center.size(); vi++)
                {
                    size_t index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                    
                    f << center[index] << "\t";
                }
                
                if ( center.size() < 3 )
                {
                    f << "0.0";
                }
                    
                f << endl;
            }
            
            // Write point data out. Only take the first value after coordinates
            
            f << endl;
            f << "POINT_DATA " << track->size() << endl;
            f << "SCALARS number_of_points FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;

            for ( typename Tracking<T>::track_t::iterator ti=track->begin(); ti!=track->end(); ti++ )
            {
                f << ti->points.size() << endl;
            }
            
            f.close();
        }
    }
    
    template <class T>
    void
    VisitUtils<T>::write_cluster_weight_response_vtk(const string &base_name,
                                                     const typename Cluster<T>::list &clusters,
                                                     WeightFunction<T> *w,
                                                     bool useMode)
    {
        string basename(base_name);
        
        for ( size_t ci = 0; ci < clusters.size(); ci++ )
        {
            Cluster<T> *cluster = clusters[ci];
            
            vector<T> mode = cluster->mode;
            
            // Write cluster out
            
            boost::replace_all( basename, "/", "_" );
            boost::replace_all( basename, "..", "" );
            
            string filename = basename + "_" + boost::lexical_cast<string>( cluster->id ) + ".vtk";
            
            ofstream f( filename.c_str() );
            f << fixed << setprecision(4);
            
            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Cluster weight function response" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << cluster->points.size() << " FLOAT" << endl;
            
            // Write point coordinates out as unstructured grid
            
            size_t point_dim = cluster->points[0]->coordinate.size();
            
            for ( size_t pi = 0; pi < cluster->points.size(); pi++ )
            {
                Point<T> *p = cluster->points[pi];
                
                for ( size_t vi = 0; vi < point_dim; vi++)
                {
                    size_t index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                    
                    f << p->coordinate[index] << "\t";
                }
                
                for (size_t j=p->coordinate.size(); j < 3; j++)
                {
                    f << "0.0\t";
                }
                
                f << endl;
            }
            
            // Write point data out. Only take the first value after coordinates
            
            
            T wr = useMode ? cluster->modal_weight_response(w) : cluster->average_weight_response(w);
            
            f << endl;
            f << "POINT_DATA " << cluster->points.size() << endl;
            f << "SCALARS weight FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;
            
            for ( size_t pi = 0; pi < cluster->points.size() - 1 ; pi++ )
            {
                f << wr << endl;
            }
            
            f.close();
        }
    }
    
    template <class T>
    void
    VisitUtils<T>::write_cluster_meanshift_vtk(const string &base_name,
                                               const typename Cluster<T>::list &clusters,
                                               bool use_ids,
                                               bool spatial_only)
    {
        std::string basename(base_name);
        boost::replace_all( basename, "/", "_" );
        boost::replace_all( basename, "..", "" );
        
        for ( size_t ci = 0; ci < clusters.size(); ci++ )
        {
            Cluster<T> *cluster = clusters[ci];
            
            typedef vector< vector<T> > vvector_t;
            
            vvector_t origins;
            
            vvector_t shifts;
            
            for (size_t pi=0; pi < cluster->points.size(); pi++)
            {
                typename Point<T>::ptr p = cluster->points[pi];
                
                if (spatial_only)
                {
                    origins.push_back(p->coordinate);
                    
                    vector<T> shift = p->shift;
                    
                    shifts.push_back(vector<T>( &shift[0], &shift[p->coordinate.size()]));
                }
                else
                {
                    origins.push_back(p->values);
                    shifts.push_back(p->shift);
                }
            }
            
            string filename = basename + "_" + boost::lexical_cast<string>( use_ids ? cluster->id : ci ) + ".vtk";
            
            ::cfa::utils::VisitUtils<T>::write_vectors_vtk(filename,origins,shifts,"shift");
        }

    }



    
}};

#endif
