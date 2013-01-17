#include <cf-algorithms/cf-algorithms.h>

#include <netcdf>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <vector>
#include <map>

#include <radolan/radolan.h>
#include <cf-algorithms/cf-algorithms.h>
#include <meanie3D/meanie3D.h>

using namespace std;
using namespace cfa::meanshift;
using namespace Radolan;
using namespace m3D;
using namespace netCDF;

typedef vector<int> coordinate_t;

/** Create gaussian distributed variable with mean 0.5 and variance
 * 
 */
float gaussian_random( int n )
{
    float sum;
    
    for (int i=0; i<n; i++)
    {
        sum += rand();
    }
    
    return sum / ((float)n);
}

/** @returns uniform random variable (0..1)
 */
float ranf()
{
    return ((float)rand())/((float)RAND_MAX);
}

/** Normal random variate generator, courtesy of
 * ftp://ftp.taygeta.com/pub/c/boxmuller.c
 * @param m mean value
 * @param s standard deviation
 * @return random variable with mean m, standard deviation s 
 */
float box_muller(float m, float s)	
{				        
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;
    static int initialized_rand = 0;
    
    if ( !initialized_rand )
    {
        srandomdev();

        initialized_rand = 1;
    }
    
//	if (use_last)		        /* use value from previous call */
//	{
//		y1 = y2;
//		use_last = 0;
//	}
//	else
//	{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
        
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
//		use_last = 1;
//	}
    
	return( m + y1 * s );
}

/** Creates a random coordinate vector.
 * @param dimension
 * @param mean
 * @param deviation
 * @return allocated random coordinate
 */
coordinate_t *randomPoint( size_t dim, float m, float s )
{
    coordinate_t *newPoint = new std::vector<int>(dim);
    
    for ( size_t i = 0; i < dim; i++ )
    {
        newPoint->at(i) = box_muller( m, s );
    }
    
    return newPoint;
}

/** Write out an axis with dim->num_dim() grid points, values varying
 * linear between min and max
 */
void writeAxis( NcFile &file, const NcDim &dim, float min, float max )
{
    // write x-axis information. Make 0 the center of the map. 
    
    float *data = (float *) malloc( sizeof(float) * dim.getSize() );

    for ( int i=0; i < dim.getSize(); i++ )
    {
        data[i] = min + ((max - min)/dim.getSize()) * i;
    }
    
    NcVar dimData = file.addVar( dim.getName(), NcType::nc_FLOAT, dim );
    
    dimData.putAtt("value_min", NcType::nc_FLOAT, min );
    
    dimData.putAtt("value_max", NcType::nc_FLOAT, max );
    
    dimData.putVar(data);
    
    free(data);
}

/**
 */
void writeCloud(NcFile &file,
                NcVar &variable,
                vector<NcDim> dims,
                const vector<float> &center,
                const vector<float> &mean,
                const vector<float> &deviation,
                size_t cloudSize )
{
    map<NcDim,NcVar> vars;
    
    for ( size_t i=0; i<dims.size(); i++ )
    {
        vars[ dims.at(i) ] = file.getVar( dims.at(i).getName() );
    }
    
    // writing individual points here, counts are all 1
    
    vector<size_t> counts( dims.size(),1 );
    
    // allocate a cursor
    
    vector<size_t> cursor( dims.size(),0 );
    
    // start generating random points (with random values between 0 and 1)
    
    size_t numPoints = 0;
    
    do 
    {
        // genrate a random coordinate
        
        for ( size_t d = 0; d < dims.size(); d++ )
        {
            bool valid = false;
            
            while ( !valid )
            {
                NcDim dim = dims[d];
                
                float rand = box_muller( mean.at(d), deviation.at(d) );
                
                // re-transform to grid coordinates
                
                float min,max;
                
                vars[dim].getAtt("value_min").getValues( &min );
                vars[dim].getAtt("value_max").getValues( &max );
                
                long n = (long)round( (dim.getSize()-1)*(rand - min) / ( max - min ) );
                
                if ( n >=0 && n < dim.getSize() )
                {
                    cursor[d] = n;
                    
                    valid = true;
                }
            }
        }
        
        // generate a random value
        
        float value = 1.0f;
        
        variable.putVar( cursor, &value );
        
        numPoints++;
        
    } while ( numPoints < cloudSize );
}

void writeCloudND( const char *filename, size_t cloud_size, size_t gridSize, float mean, float deviation, vector<float> *centerOffset = NULL )
{
    
}

void writeCloud2D( const char *filename, size_t cloud_size, size_t gridSize, float m, float s, vector<float> *centerOffset = NULL )
{
    NcFile file( filename, NcFile::replace );
    
    // Create gaussian variable with 3 dimensions
    vector<NcDim> dims;
    
    NcDim x = file.addDim("x",gridSize);
    writeAxis( file, x, -100.0, 100.0 );
    dims.push_back( x );

    NcDim y = file.addDim("y",gridSize);
    writeAxis( file, y, -100.0, 100.0 );
    dims.push_back( y );

    NcVar gaussian = file.addVar("gaussian", NcType::nc_FLOAT, dims );
    gaussian.putAtt( "valid_min", ncFloat, 0.0f );
    gaussian.putAtt( "valid_max", ncFloat, 1.0f );
    
    // write cloud out

    vector<float> mean;
    mean.push_back( m );
    mean.push_back( m );
    
    vector<float> deviation;
    deviation.push_back( s );
    deviation.push_back( s );
    
    vector<float> center;
    
    if ( centerOffset != NULL )
    {
        center = *centerOffset;
    }
    else
    {
        center.push_back(0.0);
        center.push_back(0.0);
    }
    
    writeCloud( file, gaussian, dims, center, mean, deviation, cloud_size );
}


void writeCloud3D( const char *filename, size_t cloud_size, size_t gridSize, float m, float s, vector<float> *centerOffset = NULL )
{
    NcFile file( filename, NcFile::replace );
    
    // Create gaussian variable with 3 dimensions
    vector<NcDim> dims;
    
    NcDim x = file.addDim("x",gridSize);
    writeAxis( file, x, -100.0, 100.0 );
    dims.push_back( x );
    
    NcDim y = file.addDim("y",gridSize);
    writeAxis( file, y, -100.0, 100.0 );
    dims.push_back( y );
    
    NcDim z = file.addDim("z",gridSize);
    writeAxis( file, z, -100.0, 100.0 );
    dims.push_back( z );
    
    NcVar gaussian = file.addVar( "gaussian", NcType::nc_FLOAT, dims );
    gaussian.putAtt( "valid_min", ncFloat, 0.0f );
    gaussian.putAtt( "valid_max", ncFloat, 1.0f );
    
    // write cloud out
    
    vector<float> mean;
    mean.push_back( m );
    mean.push_back( m );
    mean.push_back( m );
    
    vector<float> deviation;
    deviation.push_back( s );
    deviation.push_back( s );
    deviation.push_back( s );
    
    vector<float> center;
    
    if ( centerOffset != NULL )
    {
        center = *centerOffset;
    }
    else
    {
        center.push_back(0.0);
        center.push_back(0.0);
        center.push_back(0.0);
    }
    
    writeCloud( file, gaussian, dims, center, mean, deviation, cloud_size );
}



//int main(int argc, char** argv)
//{
//    writeCloud2D("gaussian2D.nc", 2000, 101, 0.0, 40.0 );
//    
//    // writeCloud3D("gaussian3D.nc", 4000, 101, 0.0, 40.0 );
//   
//    return 0;
//}







