#ifndef _M3D_DefaultWeightFunction_H_
#define _M3D_DefaultWeightFunction_H_

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
    using cfa::array::MultiArray;
    
    template <class T>
    class DefaultWeightFunction : public cfa::utils::WeightFunction<T>
    {
    protected:
        
        const FeatureSpace<T>       *m_fs;
        
        map<size_t,T>       m_min;
        map<size_t,T>       m_max;
        MultiArray<T>       *m_weight;
        
        void
        build_saliency_field(const FeatureSpace<T> *fs)
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
        DefaultWeightFunction(const FeatureSpace<T> *fs)
        : m_fs(fs)
        , m_min(fs->min())
        , m_max(fs->max())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes()))
        {
            build_saliency_field(fs);
        }
        
        /** Construct the weight function, using the given values
         * for valid_min/valid_max
         * @param featurespace
         * @param map of lower bounds
         * @param map of upper bounds
         */
        DefaultWeightFunction(FeatureSpace<T> *fs,
                              const DataStore<T> *data_store,
                              const map<size_t,T> &min,
                              const map<size_t,T> &max)
        : m_fs(fs)
        , m_min(min)
        , m_max(max)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes()))
        {
            build_saliency_field(fs);
        }
        
        ~DefaultWeightFunction()
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
            
            for (size_t var_index = 0; var_index < m_fs->value_rank(); var_index++)
            {
                bool is_valid = false;
                
                T value = p->values.at(m_fs->spatial_rank()+var_index);
                
                // value scaled to [0..1]
                
                T range = (m_max[var_index] - m_min[var_index]);
                
                T var_weight = (value - m_min[var_index]) / range;
                
                if (var_weight > 1)
                {
                    cerr << "ERROR: weight function at " << p->coordinate << " returns " << var_weight << endl;
                }

                sum += var_weight; //pow(var_weight, 1.0);
            }
            
            return sum;
        }
        
        /** unfavorable, since it performs a reverse lookup, which is a very
         * time-consuming operation. Use grid points where possible.
         */
        T operator()( const vector<T> &values ) const
        {
            typename CoordinateSystem<T>::GridPoint gp
                = m_fs->coordinate_system->newGridPoint();
            
            try
            {
                m_fs->coordinate_system->reverse_lookup(values,gp);
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
