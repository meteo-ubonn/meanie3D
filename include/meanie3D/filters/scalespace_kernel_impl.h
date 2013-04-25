#ifndef _M3D_ScaleSpaceKernel_Impl_H_
#define _M3D_ScaleSpaceKernel_Impl_H_

namespace m3D {
    
    template <typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const T& t)
        : m_t(t)
        , m_gauging_factor(1.0/sqrt(2.0*M_PI*t))
        , m_values(vector<T>(0))

    {
    };

    template <typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const T& t, const vector<T> &distances)
        : m_t(t)
        , m_gauging_factor(1.0/sqrt(2.0*M_PI*t))
        , m_values(vector<T>(distances.size(),0.0))
    {
        for (size_t i=0; i<distances.size(); i++)
        {
            m_values[i] = this->value(distances[i]);
        }
    };
    
    template <typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const ScaleSpaceKernel<T> &o)
        : m_t(o.m_t)
        , m_gauging_factor(o.m_gauging_factor)
        , m_values(o.m_values)
    {
    }
    
    template <typename T>
    ScaleSpaceKernel<T>::~ScaleSpaceKernel() {};
    
    template <typename T>
    ScaleSpaceKernel<T>
    ScaleSpaceKernel<T>::operator = ( const ScaleSpaceKernel<T>& other )
    {
        return ScaleSpaceKernel<T>(other);
    }
    
#pragma mark -
#pragma mark Accessors
    
    template <typename T>
    const vector<T>&
    ScaleSpaceKernel<T>::values()
    {
        return m_values;
    };
    
    template <typename T>
    bool
    ScaleSpaceKernel<T>::isPreSampled()
    {
        return m_values.size() > 0;
    }

    
#pragma mark -
#pragma mark Calculation

    template <typename T>
    const T
    ScaleSpaceKernel<T>::value(size_t index)
    {
        return m_values[index];
    }
    
    template <typename T>
    const T
    ScaleSpaceKernel<T>::value(T r)
    {
        return m_gauging_factor * std::exp(-(r*r)/(2*m_t));
    }
    
};

#endif
