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


#ifndef M3D_COORDINATESYSTEM_H
#define M3D_COORDINATESYSTEM_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>

#include <vector>
#include <map>
#include <netcdf>

namespace m3D { 
    
    using namespace ::std;
    using namespace ::netCDF;

    /** This represents one point f in feature space F.
     */
    template <class T>
    class CoordinateSystem {
    public:

#pragma mark -
#pragma mark Definitions

        typedef vector<T*> DimensionData;

        typedef vector<T> Coordinate;

        typedef vector<int> GridPoint;

    private:

#pragma mark -
#pragma mark Member variables

        DimensionData m_dimension_data;

        vector<NcDim> m_dimensions;

        vector<NcVar> m_dimension_variables;

        vector<std::string> m_dimension_units;

        vector<T> m_resolution;

        T m_resolution_norm;

        vector<size_t> m_dimension_sizes;

#pragma mark -
#pragma mark Util methods
    
        /** Calculates min/max and resolution, reads dimension variable
         * data and gets dimension sizes. Requires the dimensions and
         * dimension variables to be set up.
         */
        void construct();


        /** Reads the complete set of variables for the dimensions (like lat(x), x(x) etc.) and stores
         * it for future use in the given map. It's the caller's responsibility to free the map after
         * use.
         *
         * @param netcdf file pointer
         * @param vector of the dimensions used
         * @param associated dimension variables (in the same order!)
         * @return map with variables as key and T* as data.
         */
        static void read_dimension_data_map(DimensionData &dim_data_map,
                const vector<NcDim> &dimensions,
                const vector<NcVar> &dimVars);

        /** Iterates over the given map and frees up the allocated T* memory blocks.
         * Afterwards, the map will be empty.
         */
        static void clear_dimension_data_map(DimensionData &data);
        
    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** Creates a coordinate system from the given netcdf file
         * assuming that the dimensions are all of the format x(x),
         * y(y) etc.
         * 
         * @param filename
         * @param dimensions
         */
        CoordinateSystem(NcFile *file, 
                         const vector<string> &dimensions);
        
        /** Construct a coordinate system from the coordinate
         * variables in a NetCDF file.
         *
         * @param list of dimensions constituting this coordinate system
         * @param list of variables containing the axis data
         */
        CoordinateSystem(const vector<NcDim> &dimensions,
                const vector<NcVar> &dimension_variables);

        // big three

        /** Copy constructor
         */
        CoordinateSystem(const CoordinateSystem& other);

        /** Copy operator
         */
        CoordinateSystem<T> operator =(const CoordinateSystem<T> &other);

        /** Destructor */
        ~CoordinateSystem();

#pragma mark -
#pragma mark Transformations

        /** Lookup.
         * @throws std::out_of_range
         */
        void lookup(const GridPoint &gridpoint, Coordinate &coordinate) const;

        /** Reverse Lookup.
         * @throws std::out_of_range
         */
        void reverse_lookup(const Coordinate &coordinate, GridPoint &gridpoint) const;

        /** Lookup.
         * Allocates a new instance of Coordinate.
         * @param grid point
         * @return coordinate
         */
        Coordinate *coordinate(GridPoint &gridpoint);

        /** Reverse lookup.
         * Allocates a new instance of GridPoint.
         * @param coordinate
         * @return grid point
         */
        GridPoint *gridpoint(Coordinate &coordinate);

#pragma mark -
#pragma mark Rounding to grid

        /** Rounds the given coordinate to the closest grid
         * coordinate
         * 
         * @param point
         * @return rounded point
         */
        vector<T>
        round_to_grid(const vector<T> &v) const;

        /** Returns the grid index of the closest grid point
         * 
         * @param point
         * @return grid index
         */
        vector<int>
        rounded_gridpoint(const vector<T> &v) const;

        /** Uses resolution to convert the given vector
         * into number of grid index points. 
         * 
         * @param vector in coordinate space
         * @return vector in grid index space
         */
        vector<int>
        to_gridpoints(const vector<T> &v) const;

#pragma mark -
#pragma mark Factory methods

        /** Returns a fresh instance of a grid point, initialized to
         * the correct length and all components are zero
         */
        GridPoint newGridPoint() const;

        /** Returns a fresh instance of a coordinate, initialized to
         * the correct length and all components are zero
         */
        Coordinate newCoordinate() const;

#pragma mark -
#pragma mark Other

        /** Returns the number of dimensions
         */
        size_t rank() const;

        /** Read-only accessor
         */
        const vector<NcDim> &dimensions() const;

        /** @return the dimension variables
         */
        const vector<NcVar> &dimension_variables() const;

        /** Get the variable holding the dimension's data
         * @param dimension
         * @return variable
         */
        NcVar dimension_variable(const NcDim &dimension) const;

        /** @returs a vector with the dimension sizes
         */
        const vector<size_t> get_dimension_sizes() const;

        /** @return vector containing the size of the grid's pixels in
         * each dimension.
         */
        const vector<T>& resolution() const;

        /** The grid resolution vector is somewhat unwieldy in some estimates. This
         * method provides an easy way of getting a crude average grid resolution.
         * @return L2 norm of the grid resolution vector.
         */
        T resolution_norm() const;

        /** @return const pointer to dimension data or <code>NULL</code> 
         * if the variable does not match any dimension variables.
         * @param dimension variable
         */
        const T*
        get_dimension_data_ptr(NcVar var) const;

        /** @return const pointer to dimension data.
         * @param dimension index
         */
        const T*
        get_dimension_data_ptr(int index) const;

        /** Transforms the given coordinate such, that all entries are in meters.
         * Obviously this only applies if the dimension variables are convertible
         * into meters. If they are not in cm, m, km, inch, foot or mile, they are
         * unchanged by this transformation.
         * @param input coordinate
         * @param output coordinate
         */
        Coordinate
        to_meters(const Coordinate &coord) const;
    };
}

#endif
