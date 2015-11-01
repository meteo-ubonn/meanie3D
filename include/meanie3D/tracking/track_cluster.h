/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef M3D_TRACK_CLUSTER_H
#define	M3D_TRACK_CLUSTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/clustering.h>

#include <fstream>
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

        // The file containing the actual data
        std::string m_filename;

        // name of the cluster's variable in external memory
        std::string m_variable_name;

        // Flag indicating if points need writing off
        // or whether they can remain in memory]
        bool m_writes_points_to_disk;

        // Flag indicating if the data is available in system
        // memory or if it needs reading from external memory
        bool m_needs_reading;

        // Buffers cluster size
        size_t m_size;

        std::vector<T> m_geometrical_center;

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
            std::ofstream f(m_filename.c_str(), ios::out);
            if (!f.is_open()) {
                cerr << "FATAL:failed to open file " << m_filename << " for writing." << endl;
                exit(EXIT_FAILURE);
            }

            // three lines per point: values, gridpoint, coordinate

            typename Point<T>::list::iterator pi;
            for (pi = points.begin(); pi != points.end(); pi++) {
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

            ifstream f(m_filename.c_str(), ios::in);
            if (!f.is_open()) {
                cerr << "FATAL:failed to open file " << m_filename << " for reading." << endl;
                exit(EXIT_FAILURE);
            }

            std::string line;
            while (getline(f, line)) {
                // values
                typename Point<T>::ptr p = PointFactory<T>::get_instance()->create();
                p->values = from_string<T>(line);

                // grid point
                if (!getline(f, line)) {
                    cerr << "FATAL:failed to read line from " << m_filename << endl;
                    exit(EXIT_FAILURE);
                }
                p->gridpoint = from_string<int>(line);

                // reconstruct coordinate from spatial range
                p->coordinate = std::vector<T>(p->values.begin(),
                        p->values.begin() + p->gridpoint.size());

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
        TrackCluster(typename Cluster<T>::ptr cluster,
                bool write_points_to_disk = false)
        : m_writes_points_to_disk(write_points_to_disk)
        , m_needs_reading(true)
        {
            this->id = cluster->id;
            this->uuid = cluster->uuid;
            this->mode = cluster->mode;
            this->m_rank = cluster->rank();
            this->m_spatial_rank = cluster->spatial_rank();
            this->m_size = cluster->size();
            this->m_geometrical_center = cluster->geometrical_center();
            std::string num_postfix = boost::lexical_cast<string>(cluster->uuid);
            // name of the external 'memory'
            m_filename = "/tmp/cluster_" + num_postfix + ".txt";
            // write out to external memory
            if (m_writes_points_to_disk) {
                this->write_points(cluster->get_points());
            }
        }

        ~TrackCluster() {
            this->clear(true);
        }

#pragma mark -
#pragma mark Public member functions

        virtual void clear(bool deletion_flag = false) {
            Cluster<T>::clear(deletion_flag);
            this->m_needs_reading = true;
        }

        typename Point<T>::list &get_points() {
            if (m_writes_points_to_disk && this->m_needs_reading) {
                this->read_points();
                this->m_needs_reading = false;
            }
            return Cluster<T>::get_points();
        }

        vector<T> geometrical_center() {
            return this->m_geometrical_center;
        }

        size_t size() {
            return this->m_size;
        }
    };
};

#endif

