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
        ScalarIndex<T>      m_weight;
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
        
        OASEWeightFunction(FeatureSpace<T> *fs)
        : m_vars(fs->variables())
        , m_weight(ScalarIndex<T>(fs->coordinate_system))
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
        
        /** Smoothing distorts the value range. In order to keep 
         * the weight calculation straight, allow users to update
         * the limits
         * @param index 
         * @param min
         * @param max
         */
        void
        update_limits(const map<size_t,T> &min, const map<size_t,T> &max)
        {
            m_min = min;
            m_max = max;
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
                
                T var_weight = 1.0;
                
                if (var.getName() == "cband_radolan_rx")
                {
                    sum += (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                }
                else if (var.getName() == "msevi_l15_ir_108")
                {
                    if (value < 30)
                    {
                        sum += (30 - value) / (30 - m_min.at(var_index));
                    }
                }
                else if (var.getName() == "linet_oase_tl")
                {
                    if (value > 0)
                    {
                        sum += (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    }
                }
                else
                {
                    // value scaled to [0..1]
                    
                    T value = p->values[p->coordinate.size()+var_index];
                    
                    var_weight = (value - m_min.at(var_index)) / ( m_max.at(var_index) - m_min.at(var_index) );
                    
                    // value^2
                    //T var_weight = ( sample->at(index)->values[weight_var_index] ) * ( sample->at(index)->values[weight_var_index] );
                    
                    // value^3
                    //T value = sample->at(index)->values[weight_var_index];
                    //T var_weight = pow( value, 3 );
                    
                    // value^t
                    //T var_weight = pow( sample->at(index)->values[weight_var_index], this->feature_space->scale() );
                    
                    // 10^(value/10)
                    // T var_weight = pow( 10, values[coordinate.size()] / 10 );
                }
            }
            
            return sum;
        }
        
        /** unfavorable, since it performs a reverse lookup, which is a very
         * time-consuming operation. Use grid points where possible.
         */
        T operator()( const typename CoordinateSystem<T>::Coordinate &coordinate, const vector<T> &values ) const
        {
            return m_weight.get(coordinate);
        }
        
        T operator()( const typename CoordinateSystem<T>::GridPoint &gp, const vector<T> &values ) const
        {
            return m_weight.get(gp);
        }
        
    };

    
}};

#endif
