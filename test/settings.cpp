//
//  fs_testsettings.cpp
//  cf-algorithms
//
//  Created by Jürgen Simon on 04.05.12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include "settings.h"

FSTestSettings::FSTestSettings(size_t num_dims,
        size_t num_vars,
        vector<float> axis_bound_values,
        size_t num_gridpoints,
        std::string test_filename)
: m_num_dimensions(num_dims),
m_num_vars(num_vars),
m_num_gridpoints(num_gridpoints),
m_axis_bound_values(axis_bound_values),
m_test_filename(test_filename) { }

FSTestSettings::FSTestSettings(size_t num_dims,
        size_t num_vars,
        size_t num_gridpoints,
        std::string test_filename)
: m_num_dimensions(num_dims),
m_num_vars(num_vars),
m_num_gridpoints(num_gridpoints),
m_test_filename(test_filename)
{
    for (size_t i = 0; i < num_dims + num_vars; i++)
    {
        m_axis_bound_values.push_back(1.0);
    }
}

size_t FSTestSettings::num_dimensions()
{
    return m_num_dimensions;
}

size_t FSTestSettings::num_vars()
{
    return m_num_vars;
}

vector<float>& FSTestSettings::axis_bound_values()
{
    return m_axis_bound_values;
}

void FSTestSettings::set_axis_bound_values(vector<float> &newValue)
{
    m_axis_bound_values = newValue; // copy
}

size_t FSTestSettings::num_gridpoints()
{
    return m_num_gridpoints;
}

size_t FSTestSettings::fs_dim()
{
    return num_dimensions() + num_vars();
}

std::string FSTestSettings::test_filename()
{
    return m_test_filename;
}

void FSTestSettings::set_test_filename(std::string value)
{
    m_test_filename = value;
}
