#ifndef M3D_VARIABLEWEIGHED_H
#define M3D_VARIABLEWEIGHED_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <map>

namespace m3D {
    
    /** Function that simply creates a weight by returning
     * the value of a variable at each point. 
     * @deprecated to be removed shortly
     */
    template <class T>
    class VariableWeighed : public WeightFunction<T> {
    private:

        NcVar m_variable;

        T m_min;

        T m_range;

        void *m_data;

        vector<size_t> m_cursor;

        CoordinateSystem<T> *m_coordinate_system;

    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** Default constructor.
         * @param file NetCDF file. Must comply to cf-metadata convention.
         * @param coordinate system
         * @param variable NetCDF variable to use as weight
         */
        VariableWeighed(NcFile* file,
                CoordinateSystem<T> *coord_system,
                const NcVar &variable);

        /** Destructor */
        ~VariableWeighed();

#pragma mark -
#pragma mark Weight Function

        T operator()(const vector<T> &values) const;

        T operator()(const typename Point<T>::ptr p) const;

        T operator()(const vector<int> &gridpoint) const;

    };
}

#endif
