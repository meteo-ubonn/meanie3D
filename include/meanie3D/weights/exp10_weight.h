#ifndef M3D_EXP10WEIGHTFUNCTION_H
#define M3D_EXP10WEIGHTFUNCTION_H

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
    class EXP10WeightFunction : public WeightFunction<T>
    {
    private:

        vector<string>      m_vars;     // variables for weighting
        MultiArray<T>       *m_weight;
        const CoordinateSystem<T> *m_coordinate_system;

        void
        calculate_weight_function(FeatureSpace<T> *fs)
        {
            for (size_t i=0; i < fs->points.size(); i++)
            {
                Point<T> *p = fs->points[i];

                T saliency = this->compute_weight(p);

                m_weight->set(p->gridpoint, saliency);
            }
        };

    public:

        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        EXP10WeightFunction(FeatureSpace<T> *fs,
                            const NetCDFDataStore<T> *data_store)
        : m_vars(data_store->variable_names())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(),0.0))
        , m_coordinate_system(fs->coordinate_system)
        {
            calculate_weight_function(fs);
        }

        ~EXP10WeightFunction()
        {
            if (m_weight!=NULL)
            {
                delete m_weight;
                m_weight = NULL;
            }
        }

        /** Actual weight computation happens here
         */
        T compute_weight(Point<T> *p)
        {
            T sum = 0.0;

            size_t num_vars = p->values.size() - p->coordinate.size();

            for (size_t var_index = 0; var_index < num_vars; var_index++)
            {
                T value = p->values[p->coordinate.size()+var_index];

                sum += pow(boost::numeric_cast<double>(value),10.0);
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
