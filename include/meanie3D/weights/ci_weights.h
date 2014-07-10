#ifndef _M3D_OASE_CI_WeightFunction_H_
#define _M3D_OASE_CI_WeightFunction_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <netcdf>
#include <vector>
#include <map>
#include <cf-algorithms/cf-algorithms.h>

namespace m3D { namespace weights {
    
    using namespace netCDF;
    using namespace ::m3D;
    using std::vector;
    using std::map;
    using cfa::utils::CoordinateSystem;
    
    //
    // Constants
    //
    
    // Variables used for scoring scheme
    
    static const size_t CI_WEIGHT_NUM_VARS = 5;
    static const char * CI_WEIGHT_VARS[] = {
        "msevi_l15_ir_108",
        "msevi_l15_wv_062",
        "msevi_l15_ir_134",
        "cband_radolan_rx",
        "linet_oase_tl"
    };
    
    // Shorthands used to access variables in order
    
    static const int msevi_l15_ir_108 = 0;
    static const int msevi_l15_wv_062 = 1;
    static const int msevi_l15_ir_134 = 2;
    static const int cband_radolan_rx = 3;
    static const int linet_oase_tl = 4;
    
    // Used in converting brightness temperature
    // to radiance and vica versa
    
    static const double NU = 931.7; // cm^-1
    static const double NU3 = NU * NU * NU;

    
    /** This class represents a weight function loosely based on the CI score
     * by Walker, MacKenzie Mecicalski, Jewett (2012). Only the static score
     * criteria are used (no time differences).
     *
     * It adds score by looking for radar and lighting signatures in a 5km
     * radius around a given point.
     *
     */
    template <class T>
    class OASECIWeightFunction : public cfa::utils::WeightFunction<T>
    {
        
    private:

        //
        // Members
        //
        
        NetCDFDataStore<T>          *m_data_store;
        const CoordinateSystem<T>   *m_coordinate_system;
        MultiArray<T>               *m_weight;
        vector<NcVar>               m_vars;
        
        //
        // Attributes for calculating brightness
        // temperature from the spectral radiances
        //
        
        map<size_t,T>       m_c1;
        map<size_t,T>       m_c2;
        map<size_t,T>       m_alpha;
        map<size_t,T>       m_beta;
        
        //
        // Members for range based weight calculations
        //
        
        PointIndex<T>       *m_index;           // index for range search
        SearchParameters    *m_search_params;   // search params for search

    public:
        
        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        OASECIWeightFunction(FeatureSpace<T> *fs,
                             const std::string &filename,
                             const int time_index = -1)
        : m_data_store(NULL)
        , m_coordinate_system(fs->coordinate_system)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        {
            // Obtain a data store with all the relevant variables
            
            try
            {
                NcFile file(filename.c_str(), NcFile::read);
            
                for (size_t i=0; i < CI_WEIGHT_NUM_VARS; i++)
                {
                    NcVar var = file.getVar(CI_WEIGHT_VARS[i]);
                    
                    if (var.isNull())
                    {
                        cerr << "CRITICAL: file requires variable " << CI_WEIGHT_VARS[i] << " for CI interest weight" << endl;
                        exit(-1);
                    }
                    
                    this->m_vars.push_back(var);
                    
                    // Obtain the constants for transforming radiances
                    // into brightness temperatures fromt he attributes:
                    
                    namespace nu = cfa::utils::netcdf;
                    
                    switch (i)
                    {
                        case msevi_l15_ir_108:
                        case msevi_l15_wv_062:
                        case msevi_l15_ir_134:
                        {
                            m_c1[i] = nu::get_attribute_value<T>(var,"rad_const1");
                            m_c2[i] = nu::get_attribute_value<T>(var,"rad_const2");
                            m_alpha[i] = nu::get_attribute_value<T>(var,"alpha");
                            m_beta[i] = nu::get_attribute_value<T>(var,"beta");
                        }
                    }
                    
                }
            }
            catch (netCDF::exceptions::NcException &e)
            {
                cerr << "CRITICAL: can not read from netcdf file " << filename << endl;
                exit(-1);
            }

            // Create the data store
            
            this->m_data_store = new NetCDFDataStore<T>(filename,fs->coordinate_system,this->m_vars,time_index);
            
            // calculate the entire function as one
            
            calculate_weight_function(fs);
        }
        
        ~OASECIWeightFunction()
        {
            if (this->m_data_store != NULL)
            {
                delete this->m_data_store;
                this->m_data_store = NULL;
            }
            
            if (this->m_index != NULL)
            {
                delete this->m_index;
                this->m_index = NULL;
            }
            
            if (this->m_search_params != NULL)
            {
                delete this->m_search_params;
                this->m_search_params = NULL;
            }
        }
        
    private:
        
        void
        calculate_weight_function(FeatureSpace<T> *fs)
        {
            // Create a purely spatial index
            
            size_t spatial_dim = fs->coordinate_system->rank();
            
            vector<size_t> indexes(spatial_dim);
            
            for (size_t i=0; i<spatial_dim; i++) {
                indexes[i]=i;
            }
            
            m_index = PointIndex<T>::create( &fs->points, indexes );
            
            // 5km search radius
            
            vector<T> bandwidth(spatial_dim,5.0);
            
            m_search_params = new RangeSearchParams<T>(bandwidth);
            
            // compute the weights
            
            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];
                
                T saliency = this->compute_weight(p);
                
                m_weight->set(p->gridpoint, saliency);
            }
            
            // dispose of stuff we do not longer need
            
            delete this->m_index;
            this->m_index = NULL;
            
            delete this->m_search_params;
            this->m_search_params = NULL;
            
            delete this->m_data_store;
            this->m_data_store = NULL;
        };

    public:
        
        /** Calculates the brightness temperature in degree centigrade
         * from the given seviri count
         * @param one of msevi_l15_ir_108, msevi_l15_wv_062, msevi_l15_ir_134
         * @param radiance value for the given channel
         * @return brightness temperature in [C]
         */
        T brightness_temperature(const size_t var_index, const T &radiance)
        {
            T Tbb = m_c2[var_index] * NU / log(1 + NU3 * m_c1[var_index] / radiance );
            
            T Tb = (Tbb - m_beta[var_index]) / m_alpha[var_index];
            
            return Tb - 273.15;
        }
        
        /** Inverse calculation. Spectral radiance from temperature
         * in degree centigrade
         * @param one of msevi_l15_ir_108, msevi_l15_wv_062, msevi_l15_ir_134
         * @param brightness temperature in [C]
         * @return radiance value for the given channel
         */
        T spectral_radiance(const size_t var_index, const T &temperature)
        {
            T Tbb = (temperature + 273.15) * m_alpha[var_index] + m_beta[var_index];
            
            return NU3 * m_c1[var_index] / (exp(m_c2[var_index] * NU / Tbb) - 1);
        }
        
    private:
        
        /** Calculates the weight at the given point using the 
         * the scoring scheme.
         */
        T compute_weight(Point<T> *p)
        {
            vector<int> g = p->gridpoint;

            bool isValid = false;
            
            // TODO: validity checks
            
            T ir_108_radiance = this->m_data_store->get(msevi_l15_ir_108,g,isValid);
            T ir_108_temp = brightness_temperature(msevi_l15_ir_108,ir_108_radiance);
            
            T wv_062_rad = this->m_data_store->get(msevi_l15_wv_062,g,isValid);
            T wv_062_temp = brightness_temperature(msevi_l15_wv_062,wv_062_rad);
            
            T ir_134_rad = this->m_data_store->get(msevi_l15_ir_134,g,isValid);
            T ir_134_temp = brightness_temperature(msevi_l15_ir_134,ir_134_rad);
            
            // Calculate score
            
            int score = 0;
            
            // IR 10.7 critical value
            
            if (ir_108_temp <= 0.0)
            {
                score++;
            }
            
            // IR 0.65 - IR 10.7
            
            T delta_wv_062_ir_108 = wv_062_temp - ir_108_temp;
            
            if (delta_wv_062_ir_108 >= -35.0 && delta_wv_062_ir_108 <= -10.0)
            {
                score++;
            }
            
            // IR 13.3 - IR 10.7
            
            T delta_ir_134_ir_108 = ir_134_temp - ir_108_temp;
            
            if (delta_ir_134_ir_108 >= -25.0 && delta_ir_134_ir_108 <= -5.0)
            {
                score++;
            }
            
            // Is radar signature > 35dBZ and/or lightning present in 5km radius?
            
            bool has_lightning = false;
            bool has_radar = false;

            typename Point<T>::list *neighbors = m_index->search(p->coordinate,m_search_params);
            
            for (size_t pi = 0; pi < neighbors->size() && !(has_lightning && has_radar); pi++)
            {
                typename Point<T>::ptr n = neighbors->at(pi);
                
                bool neighbour_is_valid = false;

                if (!has_radar)
                {
                    T cband_rx = this->m_data_store->get(cband_radolan_rx,
                                                         n->gridpoint,
                                                         neighbour_is_valid);
                    
                    has_radar = (neighbour_is_valid && cband_rx >= 35.0);
                }
                
                if (!has_lightning)
                {
                    neighbour_is_valid = false;

                    T linet_count = this->m_data_store->get(linet_oase_tl,
                                                            n->gridpoint,
                                                            neighbour_is_valid);
                    
                    has_lightning = (neighbour_is_valid && linet_count > 0.0);
                }
            }
            
            delete neighbors;

            
            // If radar >= 35dBZ is present in 5km radius:
            // increase score
            
            if (has_radar)
            {
                score++;
            }
            
            // If any lightning is present in 5km radius:
            // increase score

            if (has_lightning)
            {
                score++;
            }
            
            // If both lightning and radar are present:
            // increase score once more

            if (has_radar && has_lightning)
            {
                score++;
            }
            
            return score;
        }
        
        /** unfavorable, since it performs a reverse lookup, which is a very
         * time-consuming operation. Use grid points where possible.
         */
        T operator()( const vector<T> &values ) const
        {
            typename CoordinateSystem<T>::GridPoint gp = this->m_coordinate_system->newGridPoint();
            
            try
            {
                this->m_coordinate_system->reverse_lookup(values,gp);
            }
            catch (std::out_of_range& e)
            {
                cerr << "Reverse coordinate transformation failed for coordinate=" << values << endl;
                
                return 0.0;
            }
            
            return m_weight->get(gp);
        }
        
        T operator()(const typename Point<T>::ptr p) const
        {
            return m_weight->get(p->gridpoint);
        }
        
        T operator()(const vector<int> &gridpoint) const
        {
            return m_weight->get(gridpoint);
        }
        
    };

    
}}

#endif
