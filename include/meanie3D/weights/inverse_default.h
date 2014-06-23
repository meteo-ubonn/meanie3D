#ifndef _M3D_InverseDefaultWeightFunction_H_
#define _M3D_InverseDefaultWeightFunction_H_

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
    
    template <class T>
    class InverseDefaultWeightFunction : public cfa::utils::WeightFunction<T>
    {
    protected:
        
        vector<NcVar>       m_vars;     // variables for weighting
        map<size_t,T>       m_min;      // [index,min]
        map<size_t,T>       m_max;      // [index,max]
        MultiArray<T>       *m_weight;
        CoordinateSystem<T> *m_coordinate_system;

        void
        build_saliency_field(FeatureSpace<T> *fs)
        {
            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];
                
                T saliency = compute_weight(p);
                
                m_weight->set(p->gridpoint, saliency);
            }
        };
        
    public:
        
        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        InverseDefaultWeightFunction(FeatureSpace<T> *fs,
                                     const NetCDFDataStore<T> *data_store,
                                     const map<int,double> &lower_thresholds,
                                     const map<int,double> &upper_thresholds)
        : m_vars(data_store->variables())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            // Get original limits
            
            for ( size_t index = 0; index < m_vars.size(); index++ )
            {
                T min_value,max_value;
                ::cfa::utils::netcdf::unpacked_limits<T>(m_vars[index], min_value, max_value);
                
                std::map<int,double>::const_iterator li;

                li = lower_thresholds.find(m_vars[index].getId());
                
                if (li == lower_thresholds.end())
                {
                    m_min[index] = min_value;
                }
                else
                {
                    m_min[index] = li->second;
                }
                
                li = upper_thresholds.find(m_vars[index].getId());
                
                if (li == upper_thresholds.end())
                {
                    m_max[index] = max_value;
                }
                else
                {
                    m_max[index] = li->second;
                }
            }
            
            build_saliency_field(fs);
        }
        
        /** Construct the weight function, using the given values
         * for valid_min/valid_max
         * @param featurespace
         * @param map of lower bounds
         * @param map of upper bounds
         */
        InverseDefaultWeightFunction(FeatureSpace<T> *fs,
                                     const NetCDFDataStore<T> *data_store,
                                     const map<size_t,T> &min,
                                     const map<size_t,T> &max)
        : m_vars(data_store->variables())
        , m_min(min)
        , m_max(max)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            build_saliency_field(fs);
        }
        
        ~InverseDefaultWeightFunction()
        {
            if (this->m_weight != NULL) {
                delete m_weight;
                m_weight=NULL;
            }
        }

        
        /** Actual weight computation happens here
         */
        virtual T compute_weight(Point<T> *p)
        {
            T sum = 0.0;
            
            size_t num_vars = p->values.size() - p->coordinate.size();
            
            if (p->isOriginalPoint)
            {
                for (size_t var_index = 0; var_index < num_vars; var_index++)
                {
                    NcVar var = m_vars[var_index];
                    
                    T value = p->values[p->coordinate.size()+var_index];
                    
                    // value scaled to [0..1]
                    
                    T a = - 1.0 / (m_max.at(var_index) - m_min.at(var_index));
                    
                    T b = 0.5 * (1.0 - a * (m_max.at(var_index) - m_min.at(var_index)));
                    
                    T var_weight = a * value + b;

                    sum += var_weight;
                }
            }
            
            return sum;
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
