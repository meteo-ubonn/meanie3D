#ifndef M3D_TRACK_CLUSTER_H
#define	M3D_TRACK_CLUSTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/clustering.h>

#include <exception>
#include <boost/lexical_cast.hpp>

namespace m3D {

    /** A subclass of Cluster that facilitates storing and reading
     * individual clusters from a NetCDF file. Used to make tracking
     * statistics more memory efficient.
     */
    template<typename T>
    class TrackCluster : public Cluster<T>
    {
        private:
            
#pragma mark -
#pragma mark Private members
            
            // A running number unique for each cluster 
            unsigned long       m_cid;
            
            // The file containing the actual data
            std::string         m_filename;
            
            // name of the cluster's variable in external memory
            std::string         m_variable_name;
            
            // Flag indicating if the data is available in system
            // memory or if it needs reading from external memory
            bool                m_needs_reading;
            
#pragma mark -
#pragma mark Private member functions

            /**  Write point list out to file.
             * 
             * TODO: perhaps find a faster, binary file format
             * 
             * @param points
             */
            void write_points(typename Point<T>::list points)
            {
                ofstream f(m_filename,ios::out);
                if (!f.is_open())
                {
                    cerr << "ERROR:failed to open file " << m_filename << " for writing." << endl;
                    exit(EXIT_FAILURE);
                }
                
                // three lines per point: values, gridpoint, coordinate
                
                typename Point<T>::list::iterator pi;
                for (pi=get_points().begin(); pi!=get_points().end(); pi++)
                {
                    typename Point<T>::ptr p = *pi;
                    f << p->values << endl;
                    f << p->gridpoint << endl;
                }
                
                f.close();
            }

            /** Reads points from file.
             * 
             * @param list
             */            
            void read_points()
            {
                using m3D::utils::vectors::from_string;

                this->clear(true);
                
                ifstream f(m_filename,ios::in);
                if (!f.is_open())
                {
                    cerr << "ERROR:failed to open file " << m_filename << " for reading." << endl;
                    exit(EXIT_FAILURE);
                }

                std::string line;
                while (getline(f,line))
                {
                    // values
                    typename Point<T>::ptr p = PointFactory<T>::get_instance()->create();
                    p->values = from_string<T>(line);

                    // grid point
                    if (!getline(f,line))
                    {
                        cerr << "ERROR:failed to read line from " << m_filename << endl;
                        exit(EXIT_FAILURE);
                    }
                    p->gridpoint = from_string<int>(line);

                    // reconstruct coordinate from spatial range
                    p->coordinate = std::vector<T>(p->values.begin(), 
                            p->values.begin()+p->gridpoint.size());
                    
                    this->add_point(p);
                }
            }
            
#pragma mark -
#pragma mark Constructor/Destructor

        public:

            /** Constructor
             * 
             * @param c_id this id must be unique.
             * @param cluster 
             */
            TrackCluster(unsigned long c_id, typename Cluster<T>::ptr cluster)
            : m_cid(c_id)
            , m_needs_reading(true)
            {
                this->id = cluster->id;
                this->mode = cluster->mode;
                
                std::string num_postfix = boost::lexical_cast<string>(c_id);
                
                // name of the external 'memory'
                m_filename = "/tmp/cluster_" + num_postfix + ".txt";
                
                // write out to external memory
                this->write_points(cluster->get_points());
            }
            
#pragma mark -
#pragma mark Public member functions

            void clear(bool deletion_flag=false)
            {
                Cluster<T>::clear(deletion_flag);
                this->m_needs_reading = true;
            }
            
            typename Point<T>::list &get_points()
            {
                if (this->m_needs_reading)
                {
                    this->read_points();
                    this->m_needs_reading = false;
                }
                
                return Cluster<T>::get_points();
            }
    };
};

#endif

