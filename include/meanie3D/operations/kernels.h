#ifndef M3D_KERNELS_H
#define M3D_KERNELS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

namespace m3D { 

    /** Base class for mean shift kernels. A kernel can work by either calculating it around
     * a given coordinate or - in case of isotropical, symmetrical kernels - by the distance from 
     * 0. In that case the kernel has a profile function.
     * 
     * @param T value class (float, double, etc.)
     * @param d dimensioniality of the kernel
     */

    template <class T> 
    class Kernel
    {
    protected: 

        T m_kernelSize;

    public:

        /** Default */
        Kernel() : m_kernelSize(1.0) {}

        /** If you need a specific size. 
         * Usually the bandwidth takes care of this part, though.
         */
        Kernel( T kernelSize ) : m_kernelSize( kernelSize ) {}

        virtual ~Kernel() {};

        /** Calculates the kernel value at coordinate c 
         */
        virtual T apply( const vector<T> &c ) const = 0;

        /** Calculates the kernel's <b>profile</b> function at distance d. The 
         * profile is applicable in symmetrical kernels only, where it describes
         * the kernel shape around an arbitrary point.
         */
        virtual T apply( const T dist ) const = 0;
    };

    /** Gaussian normal kernel 
     */
    template <class T> 
    class GaussianNormalKernel : public Kernel<T>
    {
    public:

        GaussianNormalKernel() : Kernel<T>() {}

        GaussianNormalKernel( T kernelSize ) : Kernel<T>(kernelSize) {}

        ~GaussianNormalKernel() {};

        T apply( const vector<T> &c ) const;

        T apply( const T dist ) const;

    };

    /** Epanechniov Kernel 
     */
    template <class T> 
    class EpanechnikovKernel : public Kernel<T>
    {
    public:

        EpanechnikovKernel() : Kernel<T>() {};

        EpanechnikovKernel( T kernelSize ) : Kernel<T>(kernelSize) {};

        ~EpanechnikovKernel() {};

        T apply( const vector<T> &c ) const;

        T apply( const T dist ) const;

    };

    /** Uniform Kernel 
     */
    template <class T> 
    class UniformKernel : public Kernel<T>
    {
    public:

        UniformKernel() : Kernel<T>() {};

        UniformKernel( T kernelSize ) : Kernel<T>( kernelSize ) {};

      T apply( const vector<T> &c ) const;

        T apply( const T dist ) const;
    };
}

#endif
