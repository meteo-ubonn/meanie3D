#ifndef M3D_OASEWEIGHTFUNCTION_H
#define M3D_OASEWEIGHTFUNCTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/weights/weight_function.h>
#include <meanie3D/index/search_parameters.h>
#include <meanie3D/operations/kernels.h>

#include <netcdf>
#include <vector>
#include <map>

namespace m3D {

    using namespace netCDF;
    using namespace ::m3D;
    using std::vector;
    using std::map;

    template <class T>
    class OASEWeightFunction : public WeightFunction<T>
    {
    private:

        vector<string>      m_vars;             // variables for weighting
        map<size_t,T>       m_min;              // [index,min]
        map<size_t,T>       m_max;              // [index,max]
        MultiArray<T>       *m_weight;
        const CoordinateSystem<T> *m_coordinate_system;

        // Weight function with range weight

        PointIndex<T>       *m_index;           // index for range search
        vector<T>           m_bandwidth;        // bandwidth for range weight
        SearchParameters    *m_search_params;   // search params for search
        Kernel<T>           *m_kernel;          // kernel for weighing 

        void
        calculate_weight_function(FeatureSpace<T> *fs)
        {
            size_t spatial_dim = fs->coordinate_system->rank();

            vector<size_t> indexes(spatial_dim);

            for (size_t i=0; i<spatial_dim; i++) indexes[i]=i;

            m_index = PointIndex<T>::create( &fs->points, indexes );

            m_kernel = new GaussianNormalKernel<T>(vector_norm(m_bandwidth));

            m_search_params = new RangeSearchParams<T>(m_bandwidth);

            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];

                T saliency = this->compute_weight(p);

                m_weight->set(p->gridpoint, saliency);
            }

            delete m_index;

            delete m_kernel;

            delete m_search_params;
        };

    public:

        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        OASEWeightFunction(FeatureSpace<T> *fs,
                           const NetCDFDataStore<T> *data_store,
                           const vector<T> &bandwidth)
        : m_vars(data_store->variable_names())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        , m_bandwidth(bandwidth)
        {
            // Get original limits

            for ( size_t index = 0; index < m_vars.size(); index++ )
            {
                m_min[index] = data_store->min(index);
                m_max[index] = data_store->max(index);
            }

            calculate_weight_function(fs);
        }

        /** Construct the weight function, using the given values
         * for valid_min/valid_max
         * @param featurespace
         * @param map of lower bounds
         * @param map of upper bounds
         */
        OASEWeightFunction(FeatureSpace<T> *fs,
                           const NetCDFDataStore<T> *data_store,
                           const vector<T> &bandwidth,
                           const map<size_t,T> &min,
                           const map<size_t,T> &max)
        : m_vars(data_store->variable_names())
        , m_min(min)
        , m_max(max)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        , m_bandwidth(bandwidth)
        {
            calculate_weight_function(fs);
        }

        ~OASEWeightFunction()
        {
            if (this->m_weight != NULL) {
                delete m_weight;
                m_weight=NULL;
            }
        }


        /** Calculate the unweighed weight at point p from 
         * cloud type, 10.8um, cband_radolan, cloud optical 
         * thickness and lightning counts
         */
        T weight_version_one(Point<T> *p)
        {
            T sum = 0.0;

            const float msevi_l2_nwcsaf_ct_multiplier = 1.0;
            const float cband_radolan_rx_multiplier = 10.0;
            const float linet_oase_tl_multiplier = 10.0;

            size_t num_vars = p->values.size() - p->coordinate.size();

            for (size_t var_index = 0; var_index < num_vars; var_index++)
            {
                string var = m_vars[var_index];

                T value = p->values[p->coordinate.size()+var_index];

                if (var == "cband_radolan_rx")
                {
                    // varies from 0 .. 1. Multiplier 10x

                    T rx_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));

                    sum += cband_radolan_rx_multiplier * rx_weight;
                }
                else if (var == "msevi_l2_cmsaf_cot")
                {
                    T cot_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += cot_weight;
                }
                else if (var == "msevi_l15_ir_108")
                {
                    T ir_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += ir_weight;
                }
                else if (var == "msevi_l2_nwcsaf_ct")
                {
//                    
//                    http://www.nwcsaf.org/HTMLContributions/CT/Prod_CT.htm
//                    0  non-processed          containing no data or corrupted data
//                    
//                    1	cloud free land 	no contamination by snow/ice covered surface, 
//                                              no contamination by clouds ; but contamination 
//                                              by thin dust/volcanic clouds not checked
//                    
//                    2	cloud free sea          no contamination by snow/ice covered surface, 
//                                              no contamination by clouds ; but contamination 
//                                              by thin dust/volcanic clouds not checked
//                    3	land contaminated by snow
//                    4	sea contaminated by snow/ice
//                    5	very low and cumuliform clouds
//                    6	very low and stratiform clouds
//                    7	low and cumuliform clouds
//                    8	low and stratiform clouds
//                    9	medium and cumuliform clouds
//                    10	medium and stratiform clouds
//                    11	high opaque and cumuliform clouds
//                    12	high opaque and stratiform clouds
//                    13	very high opaque and cumuliform clouds
//                    14	very high opaque and stratiform clouds
//                    15	high semitransparent thin clouds 	
//                    16	high semitransparent meanly thick clouds 	
//                    17	high semitransparent thick clouds 	
//                    18	high semitransparent above low or medium clouds 	
//                    19	fractional clouds (sub-pixel water clouds) 	
//                    20	undefined (undefined by CMa)
//
                    // weight: low = 1 point, medium = 1.5 points, high = 2 points, very high = 2.5 points.
                    //          stratiform = x1 cumulus = x2
                    // weight from 0 .. 5, multiplier 1x

                    float height_factor = 0.0;

                    if (value >= 7 && value <= 8)
                        height_factor = 1.0;
                    else if (value >= 9 && value <= 10)
                        height_factor = 1.5;
                    else if (value >= 11 && value <= 12)
                        height_factor = 2.0;
                    else if (value >= 13 && value <= 14)
                        height_factor = 2.5;

                    // default = stratiform
                    float type_multiplier = 1.0;

                    // double for cumuliform
                    if (value==7 || value==11 || value==13)
                    {
                        type_multiplier = 2.0;
                    }

                    T ct_weight = height_factor * type_multiplier;

                    sum += msevi_l2_nwcsaf_ct_multiplier * ct_weight;
                }
                else if (var == "linet_oase_tl")
                {
                    // varies from 0 .. 1. Multiplier 100x

                    T linet_weight = (value - m_min.at(var_index)) / (m_max.at(var_index) - m_min.at(var_index));
                    sum += linet_oase_tl_multiplier * linet_weight;
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


        T convective_initiation_weight(Point<T> *p)
        {

        }


        /** Actual weight computation happens here
         */
        T compute_weight(Point<T> *p)
        {


            typename Point<T>::list *neighbors = m_index->search(p->coordinate,m_search_params);

            T weight = 0;

            for (size_t pi = 0; pi < neighbors->size(); pi++)
            {
                typename Point<T>::ptr n = neighbors->at(pi);

                // calculate weight at that point

                T point_weight = weight_version_one(n);

                // calculate distance between n and p

                T dist = vector_norm(p->coordinate - n->coordinate);

                // calculate gaussian weighed sum 

                weight += m_kernel->apply(dist) * point_weight;
            }

            delete neighbors;

            return weight;
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
}

#endif
