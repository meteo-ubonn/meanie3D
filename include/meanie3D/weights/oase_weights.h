#ifndef _M3D_OASEWeightFunction_H_
#define _M3D_OASEWeightFunction_H_

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
    class OASEWeightFunction : public cfa::utils::WeightFunction<T>
    {
    private:
        
        vector<NcVar>       m_vars;     // variables for weighting
        map<size_t,T>       m_min;      // [index,min]
        map<size_t,T>       m_max;      // [index,max]
        cfa::utils::ScalarIndex<T,T>      m_weight;
        CoordinateSystem<T> *m_coordinate_system;
        
        void
        build_saliency_field(FeatureSpace<T> *fs)
        {
            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];
                
                T saliency = this->compute_weight(p);
                
                m_weight.set(p->gridpoint, saliency);
            }
        };
        
    public:
        
        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        OASEWeightFunction(FeatureSpace<T> *fs)
        : m_vars(fs->variables())
        , m_weight(cfa::utils::ScalarIndex<T,T>(fs->coordinate_system,0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            // Get original limits
            
            for ( size_t index = 0; index < m_vars.size(); index++ )
            {
                T min_value,max_value;
                ::cfa::utils::netcdf::unpacked_limits<T>(m_vars[index], min_value, max_value);
                m_min[index] = min_value;
                m_max[index] = max_value;
            }
            
            build_saliency_field(fs);
        }
        
        /** Construct the weight function, using the given values
         * for valid_min/valid_max
         * @param featurespace
         * @param map of lower bounds
         * @param map of upper bounds
         */
        OASEWeightFunction(FeatureSpace<T> *fs, const map<size_t,T> &min, const map<size_t,T> &max)
        : m_vars(fs->variables())
        , m_min(min)
        , m_max(max)
        , m_weight(cfa::utils::ScalarIndex<T,T>(fs->coordinate_system,0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            build_saliency_field(fs);
        }
        
        /** Actual weight computation happens here
         */
        T compute_weight(Point<T> *p)
        {
            T sum = 0.0;
            
            size_t num_vars = p->values.size() - p->coordinate.size();
            
            for (size_t var_index = 0; var_index < num_vars; var_index++)
            {
                NcVar var = m_vars[var_index];
                
                T value = p->values[p->coordinate.size()+var_index];
                
                if (var.getName() == "cband_radolan_rx")
                {
                    T rx_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += rx_weight;
                }
                else if (var.getName() == "msevi_l2_cmsaf_cot")
                {
                    T cot_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += cot_weight;
                    
                }
                else if (var.getName() == "msevi_l2_nwcsaf_ct")
                {
                    T ct_weight = (value >= 7 || value <= 14) ? 1.0 : 0.0;
                    sum += ct_weight;
                }
                else if (var.getName() == "linet_oase_tl")
                {
                    T linet_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += linet_weight;
                }
                else
                {
//                    // value scaled to [0..1]
//                    T var_weight = (value - m_min.at(var_index)) / ( m_max.at(var_index) - m_min.at(var_index) );
//                    
//                    // value^2
//                    //T var_weight = ( sample->at(index)->values[weight_var_index] ) * ( sample->at(index)->values[weight_var_index] );
//                    
//                    // value^3
//                    //T value = sample->at(index)->values[weight_var_index];
//                    //T var_weight = pow( value, 3 );
//                    
//                    // value^t
//                    //T var_weight = pow( sample->at(index)->values[weight_var_index], this->feature_space->scale() );
//                    
//                    // 10^(value/10)
//                    // T var_weight = pow( 10, value / 10 );
//                    
//                    sum += var_weight;
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
                
                return 0.0;
            }
            
            return m_weight.get(gp);
        }
        
        T operator()(const typename Point<T>::ptr p) const
        {
            return m_weight.get(p->gridpoint);
        }
        
    };

    
}}

#endif
