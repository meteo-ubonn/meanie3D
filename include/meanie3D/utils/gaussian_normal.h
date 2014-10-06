#ifndef M3D_GAUSSIAN_NORMAL
#define M3D_GAUSSIAN_NORMAL

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/featurespace.h>

namespace m3D { namespace utils {

    using namespace std;

    /** Specialization of weight function as gaussian normal function
     * around a given point. Used for test purposes.
     */
    template <class T>
    class GaussianNormal : public WeightFunction<T> {
        static T FS_VALUE_MAX;

    public:

        virtual T operator()(const vector<T> &values) const {
            static float sigma_square = 0.5;

            T norm = vectors::vector_norm(values);

            T value = FS_VALUE_MAX * (1.0 / (sqrt(2 * M_PI * sigma_square))) * exp(-0.5 * norm * norm / sigma_square);

            return value;
        }

        virtual T operator()(const typename Point<T>::ptr p) const {
            return this->operator()(p->values);
        }

        virtual T operator()(const vector<int> &gridpoint) const {
            throw "not implemented";
        }

    };

    template <typename T>
    T GaussianNormal<T>::FS_VALUE_MAX = 3.0;
}}

#endif
