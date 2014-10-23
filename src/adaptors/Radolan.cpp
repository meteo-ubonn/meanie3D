
#include <meanie3D/adaptors.h>
#include <meanie3D/exceptions.h>

#include <netcdf>
#include <iostream>
#include <radolan/radolan.h>

using namespace Radolan;

#define ADD_DIMENSION_Z 0

namespace m3D { 

    /** Converts the radolan file at path into a CF-Metadata compliant NetCDF-File.
     * @param radolanPath full path to the radolan file
     * @param netcdfPath full path to the netcdf file to be created
     * @param mode NcFile::Mode for opening the netcdf file with
     * @param omitOutside @see RDReadScan
     * @return NCFile* NetCDF-Filehandler
     * @throw CFFileConversionException
     */
    netCDF::NcFile * CFConvertRadolanFile( const char* radolanPath,
                                          const char* netcdfPath,
                                          const RDDataType *threshold,
                                          netCDF::NcFile::FileMode mode,
                                          bool omitOutside)
    throw (CFFileConversionException)
    {
        if ( mode == netCDF::NcFile::read )
        {
            throw CFFileConversionException("Mode 'ReadOnly' does not make sense");
        }
        
        RDScan *scan = RDAllocateScan();
        
        int res = RDReadScan( radolanPath, scan, omitOutside );
            
        switch ( res ) 
        {
            case -1:
                throw CFFileConversionException("Insufficient memory for reading radolan scan");
                break;
                
            case -2:
                throw CFFileConversionException("Radolan file not found");
                break;

            case -3:
                throw CFFileConversionException("I/O eror when reading radolan file");
                break;
            
            default:
                break;
        }
        
        // Radolan::RDPrintScan( scan, 10, 10 );
        
        netCDF::NcFile *file = CFConvertRadolanScan( scan, netcdfPath, threshold, mode );
        
        // CFPrintConvertedRadolanScan( file, 10, 10 );
        
        RDFreeScan( scan );
        
        return file;
    }
    
    /** Converts the radolan file at path into a CF-Metadata compliant NetCDF-File. 
     * @param radolanPath full path to the radolan file
     * @param netcdfPath full path to the netcdf file to be created
     * @param mode NcFile::Mode for opening the netcdf file with
     * @return NCFile* NetCDF-Filehandler
     * @throw CFFileConversionException
     */
    netCDF::NcFile * CFConvertRadolanScan( RDScan *scan,
                                  const char* netcdfPath,
                                  const RDDataType *threshold,
                                  netCDF::NcFile::FileMode mode)

    throw (CFFileConversionException)
    {
        using namespace netCDF;
        using namespace std;
        
        NcFile* file = NULL;
        
        try
        {
            file = new netCDF::NcFile( netcdfPath, mode );
        }
        catch (const netCDF::exceptions::NcException &e)
        {
            cerr << "Exception while creating file " << netcdfPath << " : " << e.what() << endl;
            
            throw CFFileConversionException( e.what() );
        }

        // Global attributes
        
        file->putAtt( "Conventions", "CF 1.6" );
        
        file->putAtt("title", "Radolan composite in NetCDF/CF-Metadata form." );
        
        file->putAtt("institution", "HErZ-TB1 Workgroup");
        
        file->putAtt("version", "1.0");
        
        // Dimensions
        
        vector<NcDim> dims;
        
        NcDim dimT = file->addDim( "time", 1 );

        NcDim dimX = file->addDim( "x", scan->dimLon );
        
        NcDim dimY = file->addDim( "y", scan->dimLat );

#if ADD_DIMENSION_Z
       NcDim dimZ = file->addDim( "z", 1 );
#endif
        
        //dims.push_back(dimT);
#if ADD_DIMENSION_Z
        dims.push_back(dimZ);
#endif
        dims.push_back(dimY);
        dims.push_back(dimX);

        RDCoordinateSystem rcs = RDCoordinateSystem( scan->header.scanType );

        // Coordinates
        
        netCDF::NcVar x = file->addVar( "x", ncDouble, dimX );
        
        x.putAtt( "standard_name", "projection_x_coordinate" );
        
        x.putAtt( "units", "km" );
        
        NcVar y = file->addVar( "y", ncDouble, dimY );
        
        y.putAtt( "standard_name", "projection_y_coordinate" );
        
        y.putAtt( "units", "km" );
        
#if ADD_DIMENSION_Z
        
        NcVar z = file->addVar( "z", ncDouble, dimZ );
        
        z.putAtt( "standard_name", "projection_z_coordinate" );
        
        z.putAtt( "units", "km" );
#endif
        
        // Grid Mapping
        
        RDGridPoint origin = rdGridPoint( 0, 0 );
       
        RDGeographicalPoint origin_geo = rcs.geographicalCoordinate( origin );
        
        NcVar crs = file->addVar( "crs", NcType::nc_BYTE, dims ); // note: type is of no consequence
        
        crs.putAtt( "grid_mapping_name", "polar_stereographic" );
        
        crs.putAtt( "longitude_of_projection_origin", NcType::nc_DOUBLE, origin_geo.longitude );
        
        crs.putAtt( "latitude_of_projection_origin", NcType::nc_DOUBLE, origin_geo.latitude );
        
        crs.putAtt( "false_easting", NcType::nc_DOUBLE, 0.0f );
        
        crs.putAtt( "false_northing", NcType::nc_DOUBLE, 0.0f );
        
        crs.putAtt( "scale_factor_at_projection_origin", NcType::nc_DOUBLE, rcs.polarStereographicScalingFactor( origin_geo.longitude, origin_geo.latitude) );
        
        crs.putAtt( "units", "km" );

        // Data
        
//        NcVar *data = file->add_var( "reflectivity", ncFloat, dimX, dimY, dimZ );
        
        NcVar data = file->addVar( RDScanTypeToString(scan->header.scanType), NcType::nc_FLOAT, dims );
        
        data.putAtt("units", RDUnits( scan->header.scanType ) );
        
        data.putAtt("grid_mapping", "polar_stereographic");
        
        data.putAtt("radolan_product", RDScanTypeToString( scan->header.scanType ) );
        
        data.putAtt("standard_name", CFRadolanDataStandardName( scan->header.scanType ) );
        
        // NOTE: this value is increased by 0.5 to make the lowest value (-32.5) 
        // drop out later. 
        // TODO: define a threshold parameter to feature space constructor and
        // clustering tool command line
    
        data.putAtt( "valid_min", NcType::nc_FLOAT, RDMinValue( scan->header.scanType ) );
        
        data.putAtt( "valid_max", NcType::nc_FLOAT, RDMaxValue( scan->header.scanType ) );
        
        data.putAtt( "_FillValue", NcType::nc_FLOAT, RDMissingValue(scan->header.scanType));
        
        // TIME
        
        double timestamp = (double) RDScanTimeInSecondsSinceEpoch( scan );
        
        NcVar time = file->addVar( "time", NcType::nc_DOUBLE, dimT );
        
        time.putVar( &timestamp );
        
        time.putAtt("units", "seconds since 1970-01-01 00:00:00.0" );
        
        time.putAtt("calendar", "gregorian" );

        time.putAtt("standard_name", "time" );
        
        // Re-package data
        
        // x and y are switched around in the data (following
        // the cf-metadata convention)
        
#if ADD_DIMENSION_Z
        RDDataType *converted = (RDDataType *) malloc(scan->dimLon * scan->dimLat * sizeof(RDDataType)) ;
#else
        RDDataType *converted = (RDDataType *) malloc(scan->dimLon * scan->dimLat * sizeof(RDDataType));
#endif


        // vector<size_t> index( dims.size() );
        
        for (int iy = 0; iy < scan->dimLon; iy++ )
        {
            for ( int ix = 0; ix < scan->dimLat; ix++ )
            {
                RDDataType val = scan->data[ iy * scan->dimLat + ix ];
                
                size_t converted_index = iy * scan->dimLat + ix;

                if (val == RDMissingValue(scan->header.scanType))
                {
                    converted[converted_index] = RDMissingValue(scan->header.scanType);
                }
                else
                {
                    bool should_write = (threshold==NULL) ? true : (val >= (*threshold));
                    
                    converted[converted_index] = should_write ? val : RDMinValue(scan->header.scanType);
                }
            }
        }
        
#if ADD_DIMENSION_Z
        vector<size_t> startp(3,0);
        vector<size_t> countp(3,0);

        // z
        countp[0] = 1;

        // y
        countp[1] = scan->dimLat;
        
        // x
        countp[2] = scan->dimLon;
#else
        vector<size_t> startp(2,0);
        
        vector<size_t> countp(2,0);

        // y
        countp[0] = scan->dimLat;
        
        // x
        countp[1] = scan->dimLon;
#endif
        
        // write out
        
        try
        {
            data.putVar(startp,countp,converted);
        }
        catch (const std::exception &e)
        {
            cerr << e.what() << endl;
            exit(-1);
        }
        
        // clean up
        
        delete converted;
        
        // write x-axis information
        
        float *xData = (float *) malloc( sizeof(float) * scan->dimLon );
        
        for ( int i=0; i < scan->dimLon; i++ )
        {
            RDCartesianPoint cp = rcs.cartesianCoordinate(rdGridPoint(i, 0));
            
            xData[i] = cp.x;
        }
        
        x.putVar(xData );
        
        x.putAtt( "valid_min", NcType::nc_FLOAT, xData[0] );

        x.putAtt( "valid_max", NcType::nc_FLOAT, xData[ scan->dimLon-1 ] );
        
        free(xData);
        
        
        // write y-axis information

        float *yData = (float *) malloc( sizeof(float) * scan->dimLat );
        
        for ( int i=0; i < scan->dimLat; i++ )
        {
            RDCartesianPoint cp = rcs.cartesianCoordinate(rdGridPoint(0, i));
            
            yData[i] = cp.y;
        }
        
        y.putVar( yData );
        
        y.putAtt( "valid_min", NcType::nc_FLOAT, yData[0] );
        
        y.putAtt( "valid_max", NcType::nc_FLOAT, yData[ scan->dimLon-1 ] );
        
        free(yData);
        
        // write z-Axis information
        
#if ADD_DIMENSION_Z
        float *zData = (float *) malloc( sizeof(float) * 1 );
        zData[0] = 0.0;
        z.putVar( zData );
        z.putAtt( "valid_min", NcType::nc_FLOAT, zData[0] );
        z.putAtt( "valid_max", NcType::nc_FLOAT, zData[0] );
        free( zData );
#endif
        
        return file;
    }
    
    const char * CFRadolanDataStandardName( RDScanType t )
    {
        const char *result = NULL;
        
        switch (t) 
        {
            case RD_RX:
            case RD_EX:
                result = "reflectivity";
                break;
            case RD_RZ:
            case RD_RY:
            case RD_RV:	
            case RD_EZ:
            case RD_RH:
            case RD_RB:
            case RD_RW:
            case RD_RL:
            case RD_RU:
            case RD_RS:
            case RD_RQ:
            case RD_SQ:
            case RD_SH:
            case RD_SF:
                result = "rainrate";
                break;
            default:
                result = "unknown";
                break;
        }
        return result;
    }
    
    void CFPrintConvertedRadolanScan( netCDF::NcFile *file, int latVertices, int lonVertices )
    {
        using namespace netCDF;
        
        int dimLon = file->getDim( "x" ).getSize();

        int dimLat = file->getDim( "y" ).getSize();
        
#if ADD_DIMENSION_Z
        int dimZ = file->getDim( "z" ).getSize();
#endif
        
        // int dimT = file->get_dim(NcToken("time"))->size();
        
        NcVar data = file->getVar("reflectivity");
        
        NcVarAtt product = data.getAtt("radolan_product");
        
        std::string typeIdentifier;
        
        product.getValues( typeIdentifier );
        
        RDScanType scanType = RDScanTypeFromString( typeIdentifier.c_str() );
        
#if ADD_DIMENSION_Z
        float values[dimZ][dimLat][dimLon];
#else
        float values[dimLat][dimLon];
#endif
        data.getVar( values );
        
#if ADD_DIMENSION_Z
        for ( int iz = 0; iz < dimZ; iz++ )
        {
#endif
            for ( int iy = 0; iy < dimLat; iy++ )
            {
                if ( iy % latVertices == 0 )
                {
                    for( int ix = 0; ix < dimLon; ix ++ )
                    {
                        if  ( ix % lonVertices == 0 ) 
                        {
#if ADD_DIMENSION_Z
                            float value = values[iz][iy][ix];
#else
                            float value = values[iy][ix];
#endif
                            std::cout << (RDIsCleanMeasurementAndNotMin( scanType, value ) ? "*" : " ") ;
                        }
                    }
                    std::cout << std::endl;
                }
            }
#if ADD_DIMENSION_Z
            std::cout << std::endl;
        }
#endif
    }
}
