#ifndef M3D_VariableWeighed_Impl_H
#define M3D_VariableWeighed_Impl_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils/netcdf_utils.h>

#include "variable_weighed.h"

namespace m3D {

    template <class T>
    VariableWeighed<T>::VariableWeighed(NcFile* file, CoordinateSystem<T> *coord_system, const NcVar &variable)
    : m_variable(variable), m_coordinate_system(coord_system)
    {
        T max;

        variable.getAtt("valid_min").getValues(&m_min);
        variable.getAtt("valid_max").getValues(&max);
        m_range = (max - m_min);

        // read all data at once

        T *data = (T *) calloc(utils::netcdf::num_vals(variable), sizeof (T));

        variable.getVar(data);

        m_data = (void *) data;

        m_cursor = vector<size_t>(coord_system->rank());

        for (size_t index = 0; index < coord_system->rank(); index++) {
            m_cursor[index] = 0;
        }
    }

    template <class T>
    VariableWeighed<T>::~VariableWeighed()
    {
    }

    template <class T>
    T VariableWeighed<T>::operator()(const vector<T> &values) const
    {
        return values.back();
    }

    template <class T>
    T VariableWeighed<T>::operator()(const typename Point<T>::ptr p) const
    {
        return p->values.back();
    }

    template <class T>
    T VariableWeighed<T>::operator()(const vector<int> &gridpoint) const
    {
        throw "not implemented";
    }
}

#endif