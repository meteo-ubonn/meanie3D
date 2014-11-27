#ifndef M3D_NETCDFDATASTORE_IMPL_H
#define M3D_NETCDFDATASTORE_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/array/multiarray_blitz.h>
#include <meanie3D/utils.h>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <limits>
#include <fstream>

namespace m3D { 
    
    using namespace netCDF;

    /** Specialization of DataStore for NetCDF. 
     * 
     * TODO: handle the concept of time better!
     */
    template <typename T>
    class NetCDFDataStore : public DataStore<T>
    {

    private:

        static const T NO_VALUE;

        typedef std::map< size_t, MultiArray<T> * > multiarray_map_t;

        std::string     m_filename;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is less than the given threshold,
         * will not be added.
         */
        map<int,double> m_lower_thresholds;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is more than the given threshold,
         * will not be added.
         */
        map<int,double> m_upper_thresholds;

        /** Parameter used for construction. Points in feature-space,
         * where one variable's value is out of valid range are usually
         * ignored. By providing fill values, it is possible to replace
         * the missing value with the value provided, thus allowing
         * use of that point in spite of missing values.
         */
        map<int,double> m_replacement_values;

        std::vector<std::string>    m_variable_names;

        /** Not sure this needs holding on to, but hey
         */
        const CoordinateSystem<T> *m_coordinate_system;

        /** Index of the time(time) variable to use. If -1, it is assumed
         * that there is no time variable and it is omitted.
         */
        int m_time_index;

        // Some variable attributes that are buffered for speed
        T* m_scale_factor;
        T* m_offset;
        T* m_valid_min;
        T* m_valid_max;
        T* m_fill_value;

        // "derived" boundaries
        T* m_min;
        T* m_max;

        // Map of data for construction. The data in this
        // array is already unpacked

        multiarray_map_t    m_buffered_data;

    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** TODO: either hand in variable names as strings or give NcFile pointer.
         * It makes no sense to just hand over a filename and then an array of NcVars,
         * which necessitate an open file already. Also, the NcVar has the disadvantage 
         * of tying the algorithm in with NetCDF. The variable concept should be 
         * abstracted (also the dimension concept).
         */
        NetCDFDataStore(const std::string &filename,
                        const CoordinateSystem<T> *coordinate_system,
                        const vector<std::string> &variable_names,
                        const int time_index = -1)
        : m_filename(filename)
        , m_variable_names(variable_names)
        , m_coordinate_system(coordinate_system)
        , m_time_index(time_index)
        {
            NcFile file(filename.c_str(), NcFile::read);

            if (file.isNull())
            {
                throw "ERROR:could not open file "+filename;
            }

            m_scale_factor = new T[variable_names.size()];
            m_offset = new T[variable_names.size()];
            m_valid_min = new T[variable_names.size()];
            m_valid_max = new T[variable_names.size()];
            m_min = new T[variable_names.size()];
            m_max = new T[variable_names.size()];
            m_fill_value = new T[variable_names.size()];

            // Check if the variables exist

            for (int i=0; i<variable_names.size(); i++)
            {
                try
                {
                    NcVar var = file.getVar(variable_names[i]);

                    if (var.isNull())
                    {
                        cerr << "FATAL: no variable "+variable_names[i]+" found in file " + filename << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                catch (netCDF::exceptions::NcException &e)
                {
                    cerr << "FATAL: can't access variable "+variable_names[i]+" in file " + filename << endl;
                    exit(EXIT_FAILURE);
                }

                m_scale_factor[i] = 1.0;
                m_offset[i] = 0.0;
                m_fill_value[i] = NO_VALUE;

                m_valid_min[i] = std::numeric_limits<T>::min();
                m_min[i] = std::numeric_limits<T>::min();

                m_valid_max[i] = std::numeric_limits<T>::max();
                m_max[i] = std::numeric_limits<T>::max();
            }

            this->read();
        }

        /** Destructor
         */
        ~NetCDFDataStore()
        {
            this->discard_buffer();

            delete [] m_offset;
            delete [] m_scale_factor;
            delete [] m_fill_value;
            delete [] m_valid_min;
            delete [] m_valid_max;
            delete [] m_min;
            delete [] m_max;
        }

#pragma mark -
#pragma mark Accessors

        /** @return coordinate system */
        const CoordinateSystem<T> *coordinate_system() const
        {
            return m_coordinate_system;
        }

        const std::vector<std::string> &variable_names() const {
            return m_variable_names;
        }

        /** @return valid min attributes
         */
        const T valid_min(size_t index) const { return m_valid_min[index]; }

        /** @return valid max attributes
         */
        const T valid_max(size_t index) const { return m_valid_max[index]; }

        /** @return scale factor attributes
         */
        const T scale_factor(size_t index) const { return m_scale_factor[index]; }

        /** @return add_offset attributes
         */
        const T add_offset(size_t index) const { return m_offset[index]; }

        /** @return fill_value attributes
         */
        const T fill_value(size_t index) const { return m_fill_value[index]; }

        /** @return filename 
         */
        const std::string filename() const { return m_filename; }

#pragma mark -
#pragma mark Buffered data handling

    private:

        void
        read_buffer(size_t variable_index)
        {
            bool time_is_an_issue = this->m_time_index >= 0;

            NcFile file(this->m_filename, NcFile::read);

            NcVar variable = file.getVar(m_variable_names[variable_index]);

            // If we have to take time into the picture,
            // the number of dimensions for the data is
            // one more

            int spatial_dims = time_is_an_issue
                ? variable.getDimCount() - 1
                : variable.getDimCount();

            vector<size_t> dims(spatial_dims);
            for (size_t i=0; i<spatial_dims; i++)
            {
                dims[i] = coordinate_system()->get_dimension_sizes()[i];
            }

            MultiArray<T> *index = new MultiArrayBlitz<T>(dims,0);

            vector<int> gp(spatial_dims);

            if (spatial_dims==1)
            {
                // 1D

                int N = time_is_an_issue
                    ? variable.getDim(1).getSize()
                    : variable.getDim(0).getSize();

                T *values = (T*) calloc(N,sizeof(T));

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue) start[0] = this->m_time_index;

#if DEBUG_NETCDF_DATASTORE
                cout << "1D data set: N=" << N << endl;
#endif

                // how much to read?

                vector<size_t> count(variable.getDimCount(),0);

                if (time_is_an_issue)
                {
                    count[0] = 1;
                    count[1] = N;
                }
                else
                {
                    count[0] = N;
                }

                // read the chunk

                variable.getVar(start, count, values);

                for (int i = 0; i < N; i++)
                {
                    gp[0] = i;
                    index->set(gp,values[i]);
                }

                // Clean up

                delete values;
            }
            else if (spatial_dims == 2)
            {
                // 2D

                int N,M;

                if (time_is_an_issue)
                {
                    N = variable.getDim(1).getSize();
                    M = variable.getDim(2).getSize();
                }
                else
                {
                    N = variable.getDim(0).getSize();
                    M = variable.getDim(1).getSize();
                }
                
#if DEBUG_NETCDF_DATASTORE
                cout << "2D data set: N=" << N << " M=" << M << endl;
#endif
                // Create buffer to hold data

                T *values = (T*) calloc(N*M,sizeof(T));

                if (values==NULL) {
                    cerr << "FATAL:out of memory" << endl;
                    exit(EXIT_FAILURE);
                }

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue) start[0] = this->m_time_index;

                // how much to read?

                vector<size_t> count(variable.getDimCount(),0);

                if (time_is_an_issue)
                {
                    count[0] = 1;
                    count[1] = N;
                    count[2] = M;
                }
                else
                {
                    count[0] = N;
                    count[1] = M;
                }

                // read the chunk

                variable.getVar(start, count, values);
                
                // re-package

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        gp[0] = i;
                        gp[1] = j;
                        index->set(gp,values[i*M + j]);
                    }
                }

                // Clean up

                delete[] values;
            }

            else if (spatial_dims==3)
            {
                // 3D

                int N,M,K;

                if (time_is_an_issue)
                {
                    N = variable.getDim(1).getSize();
                    M = variable.getDim(2).getSize();
                    K = variable.getDim(3).getSize();
                }
                else
                {
                    N = variable.getDim(0).getSize();
                    M = variable.getDim(1).getSize();
                    K = variable.getDim(2).getSize();
                }
                
#if DEBUG_NETCDF_DATASTORE
                cout << "3D data set: N=" << N << " M=" << M << " K=" << K << endl;
#endif

                // Read this as a succession of 2D slices
                // Allocate buffer for 2D slices

                T *values = (T*) calloc(N*M*K,sizeof(T));
                if (values == NULL)
                {
                    cerr << "FATAL:out of memory" << endl;
                    exit(EXIT_FAILURE);
                }

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue) 
                {
                    start[0] = this->m_time_index;
                    start[1] = 0;
                }
                else
                {
                    start[0] = 0;
                }

                // how much to read?

                 vector<size_t> count;

                if (time_is_an_issue)
                {
                    count.push_back(1);
                }
                count.push_back(N);
                count.push_back(M);
                count.push_back(K);

                // read the chunk

                variable.getVar(start, count, values);

                    // Re-package
                for (int i=0; i<N; i++)
                {
                    for (int j=0; j<M; j++)
                    {
                        for (int k=0; k < K; k++)
                        {
                            gp[0] = i;
                            gp[1] = j;
                            gp[2] = k;

                            index->set(gp, values[ M*K*i + K*j +k ]);
                        }
                    }
                }

                // Clean up
                delete values;
            }
            else
            {
                cerr << "FATAL: Variables with " << spatial_dims << " spatial dimensions are not currently handled" << endl;
                exit(EXIT_FAILURE);
            }

            // Store result in map

            this->m_buffered_data[variable_index] = index;
        }

        void
        write_buffered_data(size_t variable_index)
        {
            bool time_is_an_issue = this->m_time_index >= 0;

            // Open file for writing

            NcFile file(this->m_filename, NcFile::write);

            NcVar variable = file.getVar(m_variable_names[variable_index]);

            // If we have to take time into the picture,
            // the number of dimensions for the data is
            // one more

            int spatial_dims = time_is_an_issue
            ? variable.getDimCount() - 1
            : variable.getDimCount();

            vector<size_t> dims(spatial_dims);
            for (size_t i=0; i<spatial_dims; i++)
            {
                dims[i] = coordinate_system()->dimensions()[i].getSize();
            }

            MultiArray<T> *index = this->get_data(variable_index);

            vector<int> gp(spatial_dims);

            if (spatial_dims==1)
            {
                // 1D

                int N = time_is_an_issue ? variable.getDim(1).getSize() : variable.getDim(0).getSize();

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue) start[0] = this->m_time_index;

                // how much to write?

                vector<size_t> count(variable.getDimCount(),0);

                if (time_is_an_issue)
                {
                    count[0] = 1;
                    count[1] = N;
                }
                else
                {
                    count[0] = N;
                }

                // copy the data into the buffer

                T *values = (T*) calloc(N,sizeof(T));

                if (values==NULL) {
                    cerr << "FATAL:out of memory" << endl;
                    exit(EXIT_FAILURE);
                }

                for (int i = 0; i < N; i++)
                {
                    values[i] = index->get(gp);
                }

                // write the chunk

                variable.putVar(values);

                // Clean up

                delete[] values;
            }
            else if (spatial_dims == 2)
            {
                // 2D

                int N,M;

                if (time_is_an_issue)
                {
                    N = variable.getDim(1).getSize();
                    M = variable.getDim(2).getSize();
                }
                else
                {
                    N = variable.getDim(0).getSize();
                    M = variable.getDim(1).getSize();
                }

                // Create buffer to hold data

                T *values = (T*) calloc(N*M,sizeof(T));

                if (values==NULL) {
                    cerr << "FATAL:out of memory" << endl;
                    exit(EXIT_FAILURE);
                }

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue) start[0] = this->m_time_index;

                // how much to read?

                vector<size_t> count(variable.getDimCount(),0);

                if (time_is_an_issue)
                {
                    count[0] = 1;
                    count[1] = N;
                    count[2] = M;
                }
                else
                {
                    count[0] = N;
                    count[1] = M;
                }

                // re-package

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        gp[0] = i;
                        gp[1] = j;

                        values[i*M+j] = index->get(gp);
                    }
                }

                // write

                variable.putVar(&values[0]);

                // Clean up

                delete[] values;
            }

            else if (spatial_dims==3)
            {
                // 3D

                int N,M,K;

                if (time_is_an_issue)
                {
                    N = variable.getDim(1).getSize();
                    M = variable.getDim(2).getSize();
                    K = variable.getDim(3).getSize();
                }
                else
                {
                    N = variable.getDim(0).getSize();
                    M = variable.getDim(1).getSize();
                    K = variable.getDim(2).getSize();
                }

                // Read this as a succession of 2D slices
                // Allocate buffer for 2D slices

                T *values = (T*) calloc(N*M*K, sizeof(T));

                if (values == NULL)
                {
                    cerr << "FATAL:out of memory" << endl;
                    exit(EXIT_FAILURE);
                }

                // where to start reading?

                vector<size_t> start(variable.getDimCount(),0);

                if (time_is_an_issue)
                {
                    start[0] = this->m_time_index;
                    start[1] = 0;
                }
                else
                {
                    start[0] = 0;
                }

                // how much to read?

                vector<size_t> count;

                if (time_is_an_issue)
                {
                    count.push_back(1);
                }
                count.push_back(N);
                count.push_back(M);
                count.push_back(K);

                // Re-package

                for (int i=0; i<N; i++)
                {
                    for (int j=0; j<M; j++)
                    {
                        for (int k=0; k < K; k++)
                        {
                            gp[0] = i;
                            gp[1] = j;
                            gp[2] = k;

                            values[ M*K*i + K*j +k ] = index->get(gp);
                        }
                    }
                }

                // write

                variable.putVar(values);

                // Clean up

                delete[] values;
            }
            else
            {
                cerr << "FATAL: Variables with " << spatial_dims << " spatial dimensions are not currently handled" << endl;
                exit(EXIT_FAILURE);
            }
        }

    public:

        void discard_buffer()
        {
            typename multiarray_map_t::iterator i;

            for (i = m_buffered_data.begin(); i != m_buffered_data.end(); ++i)
            {
                MultiArray<T> *indexPtr = i->second;
                i->second = NULL;
                delete indexPtr;
            }
        }

        void
        get_limits(size_t var_index)
        {
            using namespace std;
            using namespace netCDF;

            NcFile file(this->m_filename, NcFile::read);

            NcVar variable = file.getVar(m_variable_names[var_index]);

            // scale_factor

            T scale_factor = 1.0;

            map< std::string, NcVarAtt > attributes = variable.getAtts();

            map< std::string, NcVarAtt >::iterator fi;

            fi = attributes.find("scale_factor");

            if ( fi != attributes.end() )
            {
                fi->second.getValues( &scale_factor );
            }

            this->m_scale_factor[var_index] = scale_factor;

            // add_offset

            T offset = 0.0;

            fi = attributes.find("add_offset");

            if ( fi != attributes.end() )
            {
                fi->second.getValues( &offset );
            }

            this->m_offset[var_index] = offset;

            // Valid min/max

            //
            // NOTE: these values of valid_min/valid_max are UNPACKED
            //

            double valid_min, valid_max;

            variable.getAtt("valid_min").getValues( &valid_min );

            this->m_valid_min[var_index] = valid_min;

            this->m_min[var_index] = valid_min;

            variable.getAtt("valid_max").getValues( &valid_max );

            this->m_valid_max[var_index] = valid_max;

            this->m_max[var_index] = valid_max;

            // _FillValue

            T fillValue = 0.0;

            fi = attributes.find("_FillValue");

            if ( fi != attributes.end() )
            {
                fi->second.getValues( &fillValue );

                this->m_fill_value[var_index] = fillValue;
            }
        }

        void
        read()
        {
            for (size_t var_index=0; var_index < this->m_variable_names.size(); var_index++)
            {
                this->get_limits(var_index);

                this->read_buffer(var_index);
            }
        }

        void
        save()
        {
            for (size_t var_index=0; var_index < this->m_variable_names.size(); var_index++)
            {
                this->write_buffered_data(var_index);
            }
        }

        /** Writes a copy of the data store to disk.
         * @param new filename
         * @throws
         */
        void
        save_as(std::string filename)
        {
            if (filename.empty())
            {
                throw std::runtime_error("no available filename");
            }

            else if (!m_filename.empty() && filename.empty())
            {
                this->save();
            }
            else if (!m_filename.empty())
            {
                // If filename changed, copy the old file to the
                // new filename to preserve everything possible and
                // then re-open the copied file and save the current
                // data

                if (m_filename != filename )
                {
                    // copy file

                    boost::filesystem::copy_file(m_filename,filename,boost::filesystem::copy_option::overwrite_if_exists);

                    this->m_filename = std::string(filename);

                    this->save();
                }
                else
                {
                    throw runtime_error("save_as called with identical filename");
                }
            }
        }

#pragma mark -
#pragma mark Accessors
        
        /** get time index
         * 
         * @return the time index this store was constructed with
         */
        int get_time_index() {
            return this->m_time_index;
        }

        /** Retrieves a value from the memory buffers
         * @param variable index
         * @param position
         * @param boolean flag, indicating whether the value 
         *        is valid by cf-metadata standards
         * @return (unpacked) value
         */
        T get(size_t variable_index, const vector<int> &gridpoint, bool &is_valid) const
        {
            T value = 0.0;
            T unpacked_value = 0.0;

            {
                typename multiarray_map_t::const_iterator i;

                {
                        i = m_buffered_data.find(variable_index);
                }

                if (i != m_buffered_data.end())
                {
                    MultiArray<T> *indexPtr = i->second;

                    value = indexPtr->get(gridpoint);
                }
                else
                {
                    cerr << "FATAL: no buffered data for variable with index " << variable_index << endl;
                    exit(EXIT_FAILURE);
                }

                if (m_fill_value[variable_index] != NO_VALUE)
                {
                    is_valid = (value != m_fill_value[variable_index]);
                }
                else
                {
                    is_valid = (value >= m_valid_min[variable_index]) 
                            && (value <= m_valid_max[variable_index]);
                }

                // scale first, then offset

                 unpacked_value = m_scale_factor[variable_index] 
                         * value + m_offset[variable_index];
            }            

            return unpacked_value;
        }

       /** Gets a point by it's linear index. The index may run
         * from 0 ... (N-1) where N is the total number of points
         * in the grid.
         */
        virtual 
        T get(size_t variable_index,
              size_t index,
              bool &is_valid) const 
        {
            return 0;
        }


        /** Values in memory buffer are 'packed'. When retrieving values
         * through get(..) they are unpacked. In order to get the packed
         * version, use this method
         * @param variable index
         * @param value
         * @return packed value
         */
        T packed_value(const size_t variable_index,
                       const T unpacked_value)
        {
           T packed_value = (unpacked_value - m_offset[variable_index]) / m_scale_factor[variable_index];

            return packed_value;
        }

        /** Allows to set a value in memory. Does NOT affect file
         * content (unless you call save() or save_as())
         * @param variable index
         * @param position
         * @param (unpacked) value to set
         */
        void set(size_t variable_index, const vector<int> &gridpoint, T value)
        {
            MultiArray<T> *indexPtr = NULL;

            typename multiarray_map_t::const_iterator i;

            i = m_buffered_data.find(variable_index);

            if (i != m_buffered_data.end())
            {
                indexPtr = i->second;
            }
            else
            {
                cerr << "FATAL: no buffered data for variable with index " << variable_index << endl;
                exit(EXIT_FAILURE);
            }

            // scale first, then offset

            T packed_value = (value - m_offset[variable_index]) / m_scale_factor[variable_index];

            indexPtr->set(gridpoint,packed_value);
        }

        const size_t rank() const
        {
            return this->m_variable_names.size();
        }

        const size_t size() const 
        {
            size_t N = 1;
            for (size_t i=0; i < this->m_coordinate_system->rank(); i++)
            {
                N *= this->m_coordinate_system->get_dimension_sizes()[i];
            }

            return N;
        }

        const vector<size_t> get_dimension_sizes() const
        {
            return this->m_coordinate_system->get_dimension_sizes();
        }

        /** @return unpacked valid_min
         */
        const T min(size_t index) const {
            return m_offset[index] + m_min[index] * m_scale_factor[index];
        }

        /** @return unpacked valid_max
         */
        const T max(size_t index) const {
            return m_offset[index] + m_max[index] * m_scale_factor[index];
        }

        /** @return reference to the multiarray used to store
         * the data in
         */
        MultiArray<T> * get_data(size_t index)
        {
            typename multiarray_map_t::const_iterator i = m_buffered_data.find(index);

            if (i != m_buffered_data.end())
            {
                return i->second;
            }
            else
            {
                cerr << "FATAL: no buffered data for variable with index " << index << endl;
                exit(EXIT_FAILURE);
            }
        }

        /** Replaces the multiarray for the given variable index.
         * The existing data set is discarded
         */
        void set_data(size_t index, MultiArray<T> *data)
        {
            typename multiarray_map_t::iterator i = m_buffered_data.find(index);

            if (i != m_buffered_data.end())
            {
                MultiArray<T> *ptr = i->second;
                m_buffered_data.erase(i);
                delete ptr;
            }

            m_buffered_data[index] = data;
        }

        void
        for_each(size_t variable_index, typename DataStore<T>::ForEachFunctor *callback)
        {
            vector<int> gridpoint(this->get_dimension_sizes().size(),0);
            this->for_each_recursive(variable_index, callback, 0, gridpoint);
        }

    private:

        void
        for_each_recursive(size_t variable_index,
                           typename DataStore<T>::ForEachFunctor *callback,
                           size_t dim_index,
                           vector<int> &gridpoint)
        {
            size_t dimSize = this->get_dimension_sizes()[dim_index];

            if (dim_index < (this->get_dimension_sizes().size()-1) )
            {
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    for_each_recursive(variable_index,callback,dim_index+1,gridpoint);
                }
            }
            else
            {
                vector<int> gIter = gridpoint;

                for (size_t i=0; i<dimSize; i++)
                {
                    gIter[dim_index] = i;

                    bool is_valid = false;

                    T value = this->get(variable_index, gIter, is_valid);

                    callback->operator()(this, variable_index, gIter, value, is_valid);
                }
            }
        }

    };

    // Static initialization of template constant

    template <typename T> 
    const T NetCDFDataStore<T>::NO_VALUE = -9999999999;
}

#endif
