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

            for ( size_t vi = 0; vi < dim_count; vi++)
            {
                size_t dim_index = vi;
                
                if ( vi < VTK_DIMENSION_INDEXES.size() )
                {
                    dim_index = VTK_DIMENSION_INDEXES.empty() ? vi : VTK_DIMENSION_INDEXES[vi];
                }

                f << mode[dim_index] << "\t";
            }

            if ( mode.size() < 3 )
            {
                f << "0.0";
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
    VisitUtils<T>::write_clusters_vtk(const string &base_name, const typename Cluster<T>::list &list, const vector<T> &bandwidths, bool use_ids)
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

            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Meanshift Clustering Result" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << list[ci]->points.size() << " FLOAT" << endl;

            // Write point coordinates out as unstructured grid

            size_t point_dim = list[ci]->points[0]->coordinate.size();

            for ( size_t pi = 0; pi < list[ci]->points.size(); pi++ )
            {
                Point<T> *p = list[ci]->points[pi];

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
            f << "POINT_DATA " << list[ci]->points.size() << endl;
            f << "SCALARS cluster FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;

            for ( size_t pi = 0; pi < list[ci]->points.size(); pi++ )
            {
                Point<T> *p = list[ci]->points[pi];

                // f << p->trajectory_length << endl;

                f << p->values[ point_dim ] << endl;
            }

            f.close();
        }
    };
    
    template <class T>
    void
    VisitUtils<T>::write_clusters_vtk( typename ClusterList<T>::ptr list )
    {
        string basename = list->source_file;
        
        for ( size_t ci = 0; ci < list->size(); ci++ )
        {
            vector<T> mode = list->clusters[ci]->mode;
            
            // Write cluster out
            
            boost::replace_all( basename, "/", "_" );
            boost::replace_all( basename, "..", "" );
            
            string filename = basename + "_cluster_" + boost::lexical_cast<string>( list->clusters[ci]->id ) + ".vtk";
            
            ofstream f( filename.c_str() );
            f << fixed << setprecision(4);
            
            // Write Header
            f << "# vtk DataFile Version 3.0" << endl;
            f << "Meanshift Clustering Result" << endl;
            f << "ASCII" << endl;
            f << "DATASET UNSTRUCTURED_GRID" << endl;
            f << "POINTS " << list->clusters[ci]->points.size() << " FLOAT" << endl;
            
            // Write point coordinates out as unstructured grid
            
            size_t point_dim = list->clusters[ci]->points[0]->coordinate.size();
            
            for ( size_t pi = 0; pi < list->clusters[ci]->points.size(); pi++ )
            {
                Point<T> *p = list->clusters[ci]->points[pi];
                
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
            f << "POINT_DATA " << list->clusters[ci]->points.size() << endl;
            f << "SCALARS cluster FLOAT" << endl;
            f << "LOOKUP_TABLE default" << endl;
            
            for ( size_t pi = 0; pi < list->clusters[ci]->points.size() - 1 ; pi++ )
            {
                Point<T> *p = list->clusters[ci]->points[pi];
                
                // f << p->trajectory_length << endl;
                
                f << p->values[ point_dim ] << endl;
            }
            
            f.close();
        }
    };
    
}};

#endif
