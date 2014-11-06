#ifndef M3D_INVERSEDEFAULTWEIGHTFUNCTION_H
#define M3D_INVERSEDEFAULTWEIGHTFUNCTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/weights/weight_function.h>

#include <netcdf>
#include <vector>
#include <map>

namespace m3D {

    using namespace netCDF;
    using namespace ::m3D;
    using std::vector;
    using std::map;

    template <class T>
    class InverseDefaultWeightFunction : public WeightFunction<T>
    {
    protected:

        vector<string>      m_vars;     // variables for weighting
        map<size_t,T>       m_min;      // [index,min]
        map<size_t,T>       m_max;      // [index,max]
        MultiArray<T>       *m_weight;
        const CoordinateSystem<T> *m_coordinate_system;

        void
        calculate_weight_function(FeatureSpace<T> *fs)
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
        : m_vars(data_store->variable_names())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
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
        InverseDefaultWeightFunction(FeatureSpace<T> *fs,
                                     const NetCDFDataStore<T> *data_store,
                                     const map<size_t,T> &min,
                                     const map<size_t,T> &max)
        : m_vars(data_store->variable_names())
        , m_min(min)
        , m_max(max)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            calculate_weight_function(fs);
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

        T operator()(const typename Point<T>::ptr p) const
        {
            return m_weight->get(p->gridpoint);
        }
    };
}

#endif
