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

#ifndef M3D_DEFAULTWEIGHTFUNCTION_H
#define M3D_DEFAULTWEIGHTFUNCTION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/weights/weight_function.h>

#include <netcdf>
#include <vector>
#include <map>

namespace m3D {

    /** This weight function scales from [0..1]. The weight is calculated
     * by iterating over all variables. For each variable a weight is 
     * generated that varies from [0..1] as the variable goes from it's
     * valid_min to valid_max (linear). The values are summed up and 
     * divided by the number of variables. 
     */
    template <class T>
    class DefaultWeightFunction : public WeightFunction<T>
    {
    private:

        const FeatureSpace<T> *m_fs;
        map<size_t, T> m_min;
        map<size_t, T> m_max;
        MultiArray<T> *m_weight;

    public:

        /** Construct the weight function, using the default values
         * for valid_min/valid_max
         * @param featurespace
         */
        DefaultWeightFunction(const FeatureSpace<T> *fs)
        : m_fs(fs)
        , m_min(fs->min())
        , m_max(fs->max())
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(), 0))
        {
            calculate_weight_function(fs);
        }

        /** Construct the weight function, using the given values
         * for valid_min/valid_max
         * @param featurespace
         * @param map of lower bounds
         * @param map of upper bounds
         */
        DefaultWeightFunction(FeatureSpace<T> *fs,
                const DataStore<T> *data_store,
                const map<size_t, T> &min,
                const map<size_t, T> &max)
        : m_fs(fs)
        , m_min(min)
        , m_max(max)
        , m_weight(new MultiArrayBlitz<T>(fs->coordinate_system->get_dimension_sizes(), 0))
        {
            calculate_weight_function(fs);
        }

        ~DefaultWeightFunction()
        {
            if (this->m_weight != NULL) {
                delete m_weight;
                m_weight = NULL;
            }
        }

    private:

        void
        calculate_weight_function(const FeatureSpace<T> *fs)
        {
            //            #if WITH_OPENMP
            //            #pragma omp parallel for 
            //            #endif
            for (size_t i = 0; i < fs->points.size(); i++) {
                Point<T> *p = fs->points[i];

                T saliency = (fs->off_limits()->get(p->gridpoint)) ? 0.0 : compute_weight(p->values);

                //                #if WITH_OPENMP
                //                #pragma omp critical 
                //                #endif
                m_weight->set(p->gridpoint, saliency);
            }
        };

        /** Actual weight computation happens here
         */

        T compute_weight(const vector<T> &values) const
        {
            T sum = 0.0;

            for (size_t var_index = 0; var_index < m_fs->value_rank(); var_index++) {
                T value = values.at(m_fs->spatial_rank() + var_index);

                // value scaled to [0..1]
                T max = m_max.find(var_index)->second;
                T min = m_min.find(var_index)->second;
                T var_weight = (value - min) / (max - min);
                sum += var_weight;
            }

            return sum / ((T) m_fs->value_rank());
        }

    public:

        /** @return pre-calculated weight
         */
        T operator()(const typename Point<T>::ptr p) const
        {
            return m_weight->get(p->gridpoint);
        }
    };
}

#endif
