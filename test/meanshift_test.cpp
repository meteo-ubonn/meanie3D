#include <cf-algorithms/cf-algorithms.h>

#include <radolan/radolan.h>
#include <cf-algorithms/cf-algorithms.h>
#include <netcdfcpp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <map>

using namespace std;

void usage() {
    cout << "meanshift_test <netcdf-file>" << endl << "Meanshift test programme" << endl;
}

template <typename T> vector<T> * random_point( vector<NcVar *> variables, vector<NcDim *> dimensions )
{
    int size = variables.size();
    
    vector<T> *coordinate = new vector<T>( size );
    
    for ( int i = 0; i < size; i++ )
    {
        NcVar *var = variables.at(i);
        
        NcDim *dim = dimensions.at(i);
        
        long rnd = random() % dim->size();
        
        coordinate->at(i) = rnd;
    }
    
    return coordinate;
}

template <typename T> vector<T> *center_point( vector<NcDim *> dimensions, const NcVar* data )
{
    int size = dimensions.size();
    
    vector<T> *coordinate = new vector<T>( size );
    
    for ( int i = 0; i < dimensions.size(); i++ )
    {
        NcDim *dim = dimensions.at(i);
        
        long value = dim->size() / 2;
        
        coordinate->at(i) = value;
    }
    
    return coordinate;
}

/* ****************************************************** */
/* NetCDF Reading Tests                                   */
/* ****************************************************** */

// Test reading netcdf slice (one row)

template <class T> void readRow( NcVar *var, long* cursor, long rowSize )
{
    long *counts = (long *) calloc( sizeof(long), var->num_dims() );
    
    memset( counts, 0, var->num_dims() );
    
    counts[ var->num_dims() - 1 ] = rowSize;
    
    T* row = (T*) calloc( sizeof(T), rowSize );
    
    var->set_cur( cursor );
    
    // var->get( row, counts );
    
    var->get( row, 1, 1, 1, rowSize );
    
    for ( long i=0; i<rowSize; i++ )
    {
        cout << "(";
        
        for ( int k=0; k < var->num_dims()-1; k++ )
        {
            cout << cursor[k] << ",";
        }
        
        cout << cursor[ var->num_dims()-1 ] + i << ") : " << row[i] << endl;
    }
    
    free( row );
    
    free( counts );
}



/** write out the trajectory to a file for visualization
 */
template <typename T> void writeTrajectoryFile2D( const char* trajectoryFilename, NcFile& file, vector<NcDim *>& dimensions, vector<vector<T> *> *trajectory )
{
    // Read axis data for those dimensions that have it
    
    map< NcVar*, float* > axisData;
    
    for ( size_t i = 0; i < dimensions.size(); i++ )
    {
        NcDim *dim = dimensions.at( i );
        
        NcVar *var = file.get_var( dim->name() );
        
        if ( var )
        {
            size_t num_vals = var->num_vals();
            
            float *data = (float *)malloc( sizeof(float) * num_vals );
            
            var->get( data, var->num_vals() );
            
            axisData[ var ] = data;
        }
    }
    
    ofstream curveFile;
    
    curveFile.open ( trajectoryFilename );
    
    curveFile << "#trajectory" << endl;

    for ( size_t i = 0; i < trajectory->size(); i++ )
    {
        vector<T> *p = trajectory->at(i);
        
        for ( size_t d = 0; d < dimensions.size(); d++ )
        {
            int n = p->at(d);
            
            float value = (float)n;
            
            NcDim *dim = dimensions.at(d);
            
            NcVar *var = file.get_var(dim->name());
            
            if ( var )
            {
                map< NcVar*,float* >::iterator iter = axisData.find(var);
            
                if ( iter != axisData.end() )
                {
                    float *axis = iter->second;
                    
                    value = axis[n];
                }
            }
                
            curveFile << setiosflags(ios::fixed) << setprecision(2) << value << " ";
        }
        
        float value = (100.0/(float)trajectory->size()) * (float)i;
        
        curveFile << setiosflags(ios::fixed) << setprecision(2) << value << endl;
    }
    
    // Free axis data
    
    for ( map<NcVar*,float*>::iterator pair = axisData.begin(); pair != axisData.end(); pair++ )
    {
        free( pair->second );
    }
    
    curveFile.close();
    
}

template <typename T> void writeTrajectoryFile3D( const char* trajectoryFilename, NcFile& file, vector<NcDim *>& dimensions, vector<vector<T> *> *trajectory )
{
    // Read axis data for those dimensions that have it
    
    map< NcVar*, float* > axisData;
    
    for ( size_t i = 0; i < dimensions.size(); i++ )
    {
        NcDim *dim = dimensions.at( i );
        
        NcVar *var = file.get_var( dim->name() );
        
        if ( var )
        {
            size_t num_vals = var->num_vals();
            
            float *data = (float *)malloc( sizeof(float) * num_vals );
            
            var->get( data, var->num_vals() );
            
            axisData[ var ] = data;
        }
    }
    
    ofstream curveFile;
    
    curveFile.open (trajectoryFilename);
    
    curveFile << "x y z trajectory" << endl;
    
    for ( size_t i = 0; i < trajectory->size(); i++ )
    {
        vector<T> *p = trajectory->at(i);
        
        for ( size_t d = 0; d < dimensions.size(); d++ )
        {
            int n = p->at(d);
            
            float value = (float)n;
            
            NcDim *dim = dimensions.at(d);
            
            NcVar *var = file.get_var(dim->name());
            
            if ( var )
            {
                map< NcVar*,float* >::iterator iter = axisData.find(var);
                
                if ( iter != axisData.end() )
                {
                    float *axis = iter->second;
                    
                    value = axis[n];
                }
            }
            
            curveFile << setiosflags(ios::fixed) << setprecision(2) << value << " ";
        }
        
        float value = (100.0/(float)trajectory->size()) * (float)i;
        
        curveFile << setiosflags(ios::fixed) << setprecision(2) << value << endl;
    }
    
    // Free axis data
    
    for ( map<NcVar*,float*>::iterator pair = axisData.begin(); pair != axisData.end(); pair++ )
    {
        free( pair->second );
    }
    
    curveFile.close();
}

template <typename T> void ms3d( const char* filename, const char* trajectoryFilename )
{
    try 
    {
        // read input file
        
        NcFile file( filename );
        
        vector<NcDim *> dimensions;
        
        vector<NcVar *> variables;
        
        // C-Band Composite, 3D+reflectivity_horizontal
        
        //        #define featureSpaceDimension 4
        //        dimensions.push_back( file.get_dim("nz") );
        //        dimensions.push_back( file.get_dim("ny") );
        //        dimensions.push_back( file.get_dim("nx") );
        //        const char* variableName = "reflectivity_horizontal";
        
        // Radolan, 2D+reflectivity
        
        //        #define featureSpaceDimension 3
        //        dimensions.push_back( file.get_dim("x") );
        //        dimensions.push_back( file.get_dim("y") );
        //        const char* variableName = "reflectivity";
        
        // 3D Gaussian Cloud
//        #define featureSpaceDimension 4
//        dimensions.push_back( file.get_dim("x") );
//        dimensions.push_back( file.get_dim("y") );
//        dimensions.push_back( file.get_dim("z") );
//        const char* variableName = "gaussian";

        
        // 3D Gaussian Cloud
        #define featureSpaceDimension 2
        dimensions.push_back( file.get_dim("x") );
        dimensions.push_back( file.get_dim("y") );
        const char* variableName = "gaussian";

        // Assemble 
        
        variables.push_back( file.get_var(variableName) );
        
        Kernel<RDDataType, featureSpaceDimension> *kernel = new GaussianNormalKernel<RDDataType, featureSpaceDimension>(100);
        
        Meanshift<RDDataType, featureSpaceDimension> ms( *kernel, dimensions, variables, 1000, 1.0/1000.0 );
        
        // run meanshift for random point with trajectory information
        
        vector<T> *randomPoint = random_point<T>( dimensions, file.get_var( variableName ) );
        
        vector<vector<T> *> *trajectory = new vector<vector<T> *>();
        
        vector<T> *destination = ms.run( &file, randomPoint, 10, trajectory );
        
        // Write to tracjectory.curve
        
//        writeTrajectoryFile3D<T>( trajectoryFilename, file, dimensions, trajectory );
        
        writeTrajectoryFile2D( trajectoryFilename, file, dimensions, trajectory );
        
//        writeTrajectoryFile3D<T>( trajectoryFilename, file, dimensions, trajectory );
        
        delete kernel;
        
        delete trajectory;
        
        delete randomPoint;
        
        delete destination;
    } 
    catch ( CFFileConversionException e )
    {
        cerr << endl << "Exception:" << e.what() << endl;
    }
}

template <typename T> void ms2d( const char* filename, const char* trajectoryFilename )
{
    try 
    {
        // read input file
        
        NcFile file( filename );
        
        vector<NcDim *> dimensions;
        
        vector<NcVar *> variables;
        
        // C-Band Composite, 3D+reflectivity_horizontal
        
        //        #define featureSpaceDimension 4
        //        dimensions.push_back( file.get_dim("nz") );
        //        dimensions.push_back( file.get_dim("ny") );
        //        dimensions.push_back( file.get_dim("nx") );
        //        const char* variableName = "reflectivity_horizontal";
        
        // Radolan, 2D+reflectivity
        
        //        #define featureSpaceDimension 3
        //        dimensions.push_back( file.get_dim("x") );
        //        dimensions.push_back( file.get_dim("y") );
        //        const char* variableName = "reflectivity";
        
        // 3D Gaussian Cloud
        //        #define featureSpaceDimension 4
        //        dimensions.push_back( file.get_dim("x") );
        //        dimensions.push_back( file.get_dim("y") );
        //        dimensions.push_back( file.get_dim("z") );
        //        const char* variableName = "gaussian";
        
        // 2D Gaussian Cloud
        #define featureSpaceDimension 2
        dimensions.push_back( file.get_dim("x") );
        dimensions.push_back( file.get_dim("y") );

        // Define feature space variables
        variables.push_back( file.get_var("x") );
        variables.push_back( file.get_var("y") );
        
        // Assemble 
        
        Kernel<RDDataType, featureSpaceDimension> *kernel = new GaussianNormalKernel<RDDataType, featureSpaceDimension>(100);
        
        Meanshift<RDDataType, featureSpaceDimension> ms( *kernel, dimensions, variables, 1000, 1.0/1000.0 );
        
        // run meanshift for random point with trajectory information
        
        vector<T> *randomPoint = random_point<T>( variables, dimensions );
        
        
        
        vector<vector<T> *> *trajectory = new vector<vector<T> *>();
        
        vector<T> *destination = ms.run( &file, randomPoint, 20, trajectory );
        
        // Write to tracjectory.curve
        
        writeTrajectoryFile2D( trajectoryFilename, file, dimensions, trajectory );
        
        //        writeTrajectoryFile3D( file, dimensions, trajectory );
        
        delete kernel;
        
        delete trajectory;
        
        delete randomPoint;
        
        delete destination;
    } 
    catch ( CFFileConversionException e )
    {
        cerr << endl << "Exception:" << e.what() << endl;
    }
}

int main(int argc, char** argv)
{
    if ( argc < 2 )
    {
        usage();
        
        exit( -1 );
    }

    // initialize random numbers
    
    srandomdev();
    
    int N = 20;

    for ( int i=0; i < N; i++ )
    {
        char fn[255];
        
        sprintf( &fn[0], "trajectory2d-%d.curve", i );
        ms2d<float>(argv[1], fn);
        
//        sprintf( &fn[0], "trajectory3d-%d.3d", i );
//        ms3d<float>( argv[1],fn );
    }

    return 0;
}