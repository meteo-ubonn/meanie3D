#ifndef M3D_COORDINATESYSTEM_IMPL_H
#define M3D_COORDINATESYSTEM_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <string.h>

#include "coordinate_system.h"

namespace m3D { 
    
    template <typename T>
    CoordinateSystem<T>::CoordinateSystem(NcFile *file, const vector<string> &dimensions)
    {
        using namespace netCDF;
        
        try
        {
            m_dimensions.clear();
            m_dimension_variables.clear();
            
            for (size_t i=0; i<dimensions.size(); i++)
            {
                m_dimensions.push_back(file->getDim(dimensions[i]));
                m_dimension_variables.push_back(file->getVar(dimensions[i]));
            }
            
            this->construct();
        }
        catch (const exceptions::NcException &e)
        {
            cerr << "ERROR:could not create coordinate system from file "
                    << file->getName()<<" :"<<e.what()<<endl;
            exit(-1);
        }
    }

    template <typename T>
    CoordinateSystem<T>::CoordinateSystem(const vector<NcDim> &dimensions,
            const vector<NcVar> &dimension_variables)
    : m_dimensions(dimensions)
    , m_dimension_variables(dimension_variables) 
    {
        this->construct();
    }

    template <typename T>
    CoordinateSystem<T>::CoordinateSystem(const CoordinateSystem& other)
    : m_dimensions(other.dimensions())
    , m_dimension_variables(other.dimension_variables())
    , m_resolution(other.m_resolution)
    , m_resolution_norm(other.m_resolution_norm) {
        // TODO: make a copy of the dimension data map, which is cheaper than reading it new

        read_dimension_data_map(m_dimension_data, m_dimensions, m_dimension_variables);

        // resolution, calculated from the dimension variables
    }

    template <typename T>
    CoordinateSystem<T>
    CoordinateSystem<T>::operator =(const CoordinateSystem<T> &other) {
        CoordinateSystem<T> copy(other);
        return copy;
    }

    template <typename T>
    CoordinateSystem<T>::~CoordinateSystem() {
        clear_dimension_data_map(m_dimension_data);
    }
    
    template <typename T>
    void
    CoordinateSystem<T>::construct()
    {
        try
        {
            // get dimension sizes

            m_dimension_sizes.resize(this->rank());

            for (size_t i = 0; i<this->rank(); i++)
                m_dimension_sizes[i] = this->m_dimensions[i].getSize();

            // Read the data for the dimensions from the dimension variables

            read_dimension_data_map(m_dimension_data, m_dimensions, m_dimension_variables);

            // resolution, calculated from the dimension variables

            m_resolution = vector<T>(m_dimension_variables.size());

            for (size_t i = 0; i < m_dimension_variables.size(); i++) 
            {
                NcVar var = m_dimension_variables[i];

                // Figure out the unit

                // TODO: expand for more units
                // TODO: test this code!!

                std::string unit_string;

                NcVarAtt unit;

                try 
                {
                    unit = var.getAtt("units");
                    unit.getValues(unit_string);
                    if (strcmp(unit_string.c_str(), "m") == 0) {
                        m_dimension_units.push_back("m");
                    } else if (strcmp(unit_string.c_str(), "km") == 0) {
                        m_dimension_units.push_back("km");
                    } else {
                        cerr << "WARNING: Variable " << var.getName() << " has units'" << unit_string << "' which is currently not handled. It will be assumed to be in meters [m]" << endl;
                        m_dimension_units.push_back("m");
                    }
                }
                catch (netCDF::exceptions::NcException &e)
                {
                    cerr << "WARNING: Variable " << var.getName() << " has no 'units' attribute. It will be assumed to be in meters [m]" << endl;
                    m_dimension_units.push_back("m");
                }

                double min, max;

                try 
                {
                    netcdf::get_valid_range(var, min, max);
                } 
                catch (std::exception &e) 
                {
                    cerr << "Variable " << var.getName() << " is missing valid_min+valid_max or valid_range attribute" << endl;
                    cerr << "Falling back on the variable values" << endl;

                    min = std::numeric_limits<T>::max();
                    max = std::numeric_limits<T>::min();

                    for (size_t vi = 0; vi < m_dimension_sizes[i]; vi++) 
                    {
                        T value = m_dimension_data[i][vi];

                        if (value > max) {
                            max = value;
                        }

                        if (value < min) {
                            min = value;
                        }
                    }

                    cerr << "Range found from data: [" << min << "," << max << "]" << endl;
                }

                m_resolution[i] = (max - min) / (m_dimension_sizes[i] - 1);
            }

            m_resolution_norm = utils::vectors::vector_norm<T>(m_resolution);
        }
        catch (const exceptions::NcException &e)
        {
            cerr << "ERROR:could not construct coordinate system:"
                    << e.what() << endl;
            exit(-1);
        }
    }

    template <typename T>
    void
    CoordinateSystem<T>::read_dimension_data_map(DimensionData &dimData,
            const vector<NcDim> &dimensions,
            const vector<NcVar> &dimVars) {
        dimData.resize(dimensions.size());

        for (size_t i = 0; i < dimVars.size(); i++) {
            NcVar var = dimVars[i];

#if DEBUG_DIMENSION_DATA
            cout << "Reading dimension data for " << var.getName() << endl;
#endif
            T* data = utils::netcdf::readNetCDFVariable<T>(var);

#if DEBUG_DIMENSION_DATA
            utils::print_array(data, utils::netcdf::num_vals(var));
            cout << endl;
#endif
            dimData[i] = data;
        }
    }

    template <typename T>
    void CoordinateSystem<T>::clear_dimension_data_map(DimensionData &data) {
        while (!data.empty()) {
            T *dataPtr = data.back();
            free(dataPtr);
            data.pop_back();
        }
    }

    template <typename T>
    typename CoordinateSystem<T>::Coordinate *
    CoordinateSystem<T>::coordinate(GridPoint &gridpoint) {
        Coordinate coordinate(gridpoint.size(), 0);

        lookup(gridpoint, coordinate);

        return new Coordinate(coordinate);
    }

    template <typename T>
    typename CoordinateSystem<T>::GridPoint *
    CoordinateSystem<T>::gridpoint(Coordinate &coordinate) {
        GridPoint gridpoint(coordinate.size(), 0);

        try {
            reverse_lookup(coordinate, gridpoint);
        } catch (std::out_of_range& e) {
            cerr << "Reverse coordinate transformation failed for coordinate=" << coordinate << endl;
        }

        return new GridPoint(gridpoint);
    }

    template <typename T>
    void CoordinateSystem<T>::lookup(const GridPoint &gridpoint, Coordinate &coordinate) const {
        assert(gridpoint.size() == coordinate.size());

        for (size_t index = 0; index < gridpoint.size(); index++)
            coordinate[index] = m_dimension_data[index][gridpoint[index]];
    }

    template <typename T>
    void CoordinateSystem<T>::reverse_lookup(const Coordinate &coordinate, GridPoint &gridpoint) const {
        vector<int> result = this->rounded_gridpoint(coordinate);

        for (size_t index = 0; index < coordinate.size(); index++) {
            if (result[index] < 0 || result[index] >= this->m_dimension_sizes[index]) {
                throw std::out_of_range("coordinate out of range");
            }

            gridpoint[index] = result[index];
        }
    }

    template <typename T>
    vector<T>
    CoordinateSystem<T>::round_to_grid(const vector<T> &v) const {
        assert(v.size() >= this->rank());

        vector<T> result = v;

        // Normalize the vector's spatial components to align with the grid

        for (size_t ci = 0; ci < m_dimensions.size(); ci++) {
#if GRID_ROUNDING_METHOD_FLOOR
            int multiplier = floor(v[ci] / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_ROUND
            int multiplier = round(v[ci] / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_CEIL
            int multiplier = ceil(v[ci] / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_RINT
            int multiplier = rint(v[ci] / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_NONE
            int multiplier = v[ci] / this->resolution()[ci];
#else
            int multiplier = ceil(v[ci] / this->resolution()[ci]);
#endif
            result[ci] = multiplier * this->resolution()[ci];
        }

        return result;
    }

    template <typename T>
    vector<int>
    CoordinateSystem<T>::rounded_gridpoint(const vector<T> &v) const {
        assert(v.size() >= this->rank());

        vector<int> result(this->rank(), 0);

        for (size_t ci = 0; ci < this->rank(); ci++) 
        {
            T corner = m_dimension_data[ci][0];

            T distance_from_corner = v[ci] - corner;

#if GRID_ROUNDING_METHOD_FLOOR
            int multiplier = floor(distance_from_corner / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_ROUND
            int multiplier = round(distance_from_corner / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_CEIL
            int multiplier = ceil(distance_from_corner / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_RINT
            int multiplier = rint(distance_from_corner / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_NONE
            int multiplier = distance_from_corner / this->resolution()[ci];
#else
            int multiplier = ceil(distance_from_corner / this->resolution()[ci]);
#endif
            result[ci] = multiplier;
        }

        return result;
    }

    template <typename T>
    vector<int>
    CoordinateSystem<T>::to_gridpoints(const vector<T> &v) const
    {
        vector<int> result(v.size());
        for (size_t ci=0; ci < this->rank(); ci++)
        {
            T value = v.at(ci);

#if GRID_ROUNDING_METHOD_FLOOR
            int number_of_points = floor(value / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_ROUND
            int number_of_points = round(value / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_CEIL
            int number_of_points = ceil(value / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_RINT
            int number_of_points = rint(value / this->resolution()[ci]);
#elif GRID_ROUNDING_METHOD_NONE
            int number_of_points = value / this->resolution()[ci];
#else
            int number_of_points = ceil(value / this->resolution()[ci]);
#endif
            result[ci] = number_of_points;
        }

        return result;
    }

    template <typename T>
    size_t CoordinateSystem<T>::rank() const {
        return m_dimensions.size();
    }

    template <typename T>
    const vector<NcDim> & CoordinateSystem<T>::dimensions() const {
        return m_dimensions;
    }

    template <typename T>
    const vector<NcVar> &
    CoordinateSystem<T>::dimension_variables() const {
        return m_dimension_variables;
    }

    template <typename T>
    NcVar
    CoordinateSystem<T>::dimension_variable(const NcDim &dim) const {
        vector<NcDim>::const_iterator it = find(m_dimensions.begin(), m_dimensions.end(), dim);

        assert(it != m_dimensions.end());

        size_t index = it - m_dimensions.begin();

        return m_dimension_variables[ index ];
    }

    template <typename T>
    const vector<T> &
    CoordinateSystem<T>::resolution() const {
        return m_resolution;
    }

    template <typename T>
    T
    CoordinateSystem<T>::resolution_norm() const {
        return m_resolution_norm;
    }

    template <typename T>
    typename CoordinateSystem<T>::GridPoint
    CoordinateSystem<T>::newGridPoint() const {
        return GridPoint(this->rank(), 0);
    }

    template <typename T>
    typename CoordinateSystem<T>::Coordinate
    CoordinateSystem<T>::newCoordinate() const {
        return Coordinate(this->rank(), 0);
    }

    template <typename T>
    const vector<size_t>
    CoordinateSystem<T>::get_dimension_sizes() const {
        return m_dimension_sizes;
    }

    template <typename T>
    const T*
    CoordinateSystem<T>::get_dimension_data_ptr(NcVar var) const {
        T* result = NULL;

        for (size_t i = 0; i < this->rank(); i++) {
            if (m_dimension_variables[i].getName() == var.getName()) {
                result = m_dimension_data[i];
                break;
            }
        }

        return result;
    }

    template <typename T>
    const T*
    CoordinateSystem<T>::get_dimension_data_ptr(int index) const {
        return m_dimension_data[index];
    }

    template <typename T>
    typename CoordinateSystem<T>::Coordinate
    CoordinateSystem<T>::to_meters(const typename CoordinateSystem<T>::Coordinate &coord) const {
        typename CoordinateSystem<T>::Coordinate transformed(this->rank(), 0);

        for (size_t i = 0; i < this->rank(); i++) {
            ::units::values::m value_in_meters;

            // TODO: expand for more units
            // TODO: test this code!!

            if (this->m_dimension_units[i] == "m") {
                value_in_meters = ::units::values::m(coord[i]);
            } else if (this->m_dimension_units[i] == "km") {
                value_in_meters = ::units::values::km(coord[i]);
            }

            transformed[i] = value_in_meters.get();
        }

        return transformed;
    }
}

#endif
