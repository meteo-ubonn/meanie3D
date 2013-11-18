//
//  meanie3D-detect
//  cf-algorithms
//
//  Created by Jürgen Lorenz Simon on 5/3/12.
//  Copyright (c) 2012 Jürgen Lorenz Simon. All rights reserved.
//

#include <boost/algorithm/string.hpp>
#include <boost/exception/info.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/program_options.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/tokenizer.hpp>

#include <meanie3D/meanie3D.h>

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <exception>
#include <locale>
#include <limits>
#include <stdlib.h>
#include <netcdf>
#include <set>

#include <shapefil.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;

/** Feature-space data type */
typedef double T;

static const double NO_SCALE = numeric_limits<double>::min();

void parse_commmandline(program_options::variables_map vm,
                        NcFile **filePtr,
                        string &source_filename,
                        string &root_directory,
                        vector<NcDim> &dimensions,
                        vector<NcVar> &dimension_variables)
{
    if ( vm.count("source") == 0 )
    {
        cerr << "Missing 'source' argument" << endl;
        
        exit( 1 );
    }
    
    source_filename = vm["source"].as<string>();
    
    try
    {
        root_directory = vm["root"].as<string>();
    }
    catch ( const boost::exception& e )
    {
        cerr << "Missing 'root' argument" << endl;
        
        exit(-1);
    }
    
    // Open NetCDF file
    
    NcFile *file = NULL;
    
    try
    {
        file = new NcFile( source_filename, NcFile::read );
    }
    catch (const netCDF::exceptions::NcException &e)
    {
        cerr << "ERROR opening file '" << source_filename << "' : " << e.what() << endl;
        exit(-1);
    }
    
    *filePtr = file;
    
    // Extract dimensions
    
    if ( vm.count("dimensions") == 0 )
    {
        cerr << "Missing parameter --dimensions" << endl;
        
        exit( 1 );
    }
    
    // parse dimension list
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    boost::char_separator<char> sep(",");
    
    tokenizer dim_tokens( vm["dimensions"].as<string>(), sep );
    
    for ( tokenizer::iterator tok_iter = dim_tokens.begin(); tok_iter != dim_tokens.end(); ++tok_iter )
    {
        const char* name = (*tok_iter).c_str();
        
        dimensions.push_back( file->getDim( name ) );
        
        dimension_variables.push_back( file->getVar( name ) );
    }
}

void add_local_topography(NcFile &mapfile)
{
    NcFile topography_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/oase-georef-1km-germany-2d-v01b.nc",NcFile::read);
    NcFile mapstuff_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/radolan_mapstuff.nc",NcFile::read);
    
    const size_t NX=510;
    const size_t NY=472;
    const size_t NZ=60;
    
    double z_fillValue = -9999.0f;

    static float topography[900][900];

    for (size_t iy=0; iy < 900; iy++)
        for (size_t ix=0; ix < 900; ix++)
            topography[iy][ix]=0.0f;
    
    topography_file.getVar("topography_srtm").getVar(&topography[0][0]);
    
    // create the output file with dims/vars
    
    nc_redef(mapfile.getId());
    
    // dims
    
    NcDim dx = mapfile.addDim("local_x", NX);
    NcDim dy = mapfile.addDim("local_y", NY);
    NcDim dz = mapfile.addDim("local_z", NZ);
    
    float x_min = -366.462f;
    float x_max = -111.962f;
    float delta_x = round( 100.0 * (x_max - x_min) / (float)NX ) / 100.0f;
    
    float y_min = -4340.645f;
    float y_max = -4105.145f;
    float delta_y = round( 100.0 * (y_max - y_min) / (float)NY ) / 100.0f;

    float z_min = 0.0f;
    float z_max = 14.75f;
    float delta_z = round( 100.0 * (z_max - z_min) / (float)NZ ) / 100.0f;
    
    // x
    
    NcVar x = mapfile.addVar("local_x", "float", "local_x");
    x.putAtt("units", "km");
    x.putAtt("long_name", "projection_x_coordinate, local composite");
    x.putAtt("valid_min", ncFloat, x_min);
    x.putAtt("valid_max", ncFloat, x_max);
    
    // y
    
    NcVar y = mapfile.addVar("local_y", "float", "local_y");
    y.putAtt("units", "km");
    y.putAtt("long_name", "projection_y_coordinate, local composite");
    y.putAtt("valid_min", ncFloat, y_min);
    y.putAtt("valid_max", ncFloat, y_max);
    
    // z
    
    NcVar z = mapfile.addVar("local_z", "float", "local_z");
    z.putAtt("units", "km");
    z.putAtt("long_name", "height above sea level, local composite");
    z.putAtt("valid_min", ncFloat, z_min);
    z.putAtt("valid_max", ncFloat, z_max);
    
    // Create 2D topo variable
    
    vector<NcDim> dims_2D;
    dims_2D.push_back(dy);
    dims_2D.push_back(dx);
    
    NcVar topo_2D = mapfile.addVar("local_topo_2D", ncDouble, dims_2D);
    topo_2D.putAtt("_FillValue", ncDouble, z_fillValue);
    topo_2D.putAtt("units", "m");
    topo_2D.putAtt("long_name", "height above sea level");
    topo_2D.putAtt("valid_min", ncDouble, -163.0);
    topo_2D.putAtt("valid_max", ncDouble, 14750.0f);
    topo_2D.putAtt("scale_factor", ncDouble, 1.0);
    topo_2D.putAtt("add_offset", ncDouble, 0.0);
    
    // Create 3D topo variable
    
    vector<NcDim> dims_3D;
    dims_3D.push_back(dz);
    dims_3D.push_back(dy);
    dims_3D.push_back(dx);
    
    NcVar topo_3D = mapfile.addVar("local_topo_3D", ncDouble, dims_3D);
    topo_3D.putAtt("_FillValue", ncDouble, z_fillValue);
    topo_3D.putAtt("units", "m");
    topo_3D.putAtt("long_name", "height above sea level");
    topo_3D.putAtt("valid_min", ncDouble, -163.0f);
    topo_3D.putAtt("valid_max", ncDouble, 14750.0f);
    topo_3D.putAtt("scale_factor", ncDouble, 1.0);
    topo_3D.putAtt("add_offset", ncDouble, 0.0);
    
    nc_enddef(mapfile.getId());

    // dim-vars
    
    float *x_data = (float *)malloc(sizeof(float)*NX);
    for (size_t i=0; i<NX; i++)
    {
        x_data[i] = x_min + i * delta_x;
    }

    float *y_data = (float *)malloc(sizeof(float)*NY);
    for (size_t i=0; i<NY; i++)
    {
        y_data[i] = y_min + i * delta_y;
    }

    float *z_data = (float *)malloc(sizeof(float)*NZ);
    for (size_t i=0; i<NZ; i++)
    {
        z_data[i] = z_min + i * delta_z;
    }
    
    try
    {
        x.putVar(&x_data[0]);
        y.putVar(&y_data[0]);
        z.putVar(&z_data[0]);
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }
    
    // Write variables
    
    static double topo_data_2D[NY][NX];
    static double topo_data_3D[NZ][NY][NX];
    
    // figure out the shift between national and local
    
    int  shift_x = (int)ceil((x_min+523.4622));
    int  shift_y = (int)ceil((y_min+4658.645f));
    
    for (size_t iy=0; iy < NY; iy++)
    {
        for (size_t ix=0; ix < NX; ix++)
        {
            // Read the z value from the topography
            // national composite is offset by dy=318 and dx=157
            // and has twice the resolution
            
            double z_m = 0.0;
            double count = 0.0;
            
            int last_iy = -1;
            int last_ix = -1;
            
            for (int ny = (shift_y + (iy-1)/2) ; ny < (shift_y + (iy+1)/2); ny++)
            {
                for (int nx = (shift_x + (ix-1)/2) ; nx < (shift_x + (ix+1)/2); nx++)
                {
                    if ((ny >=0) && (ny < 900) && (nx >=0) && (nx < 900))
                    {
                        if (ny != last_iy && nx != last_ix)
                        {
                            count = count + 1.0f;
                            z_m += topography[ny][nx];
                            last_ix = nx;
                            last_iy = ny;
                        }
                    }
                }
            }
            
            z_m = z_m / count;
            
            // Apply interpolated values
            
            topo_data_2D[iy][ix] = z_m;
            
            double z_km = z_m / 1000.0f; // [m]->[km]
            
            int z_index = (int)ceil( z_km / 0.25f );
            
            if (z_index > (NZ-1))
                z_index = (NZ-1);
            else if (z_index < 0)
                z_index = 0;
            
            for (size_t iz=0; iz<60; iz++)
            {
                topo_data_3D[iz][iy][ix] = (iz<=z_index) ? z_m : z_fillValue;
            }

//            for (size_t iz=0; iz < NZ; iz++)
//            {
//                
//                if (iz < z_index)
//                {
//                    topo_data_3D[iz][iy][ix] = z_data[iz];
//                }
//                else if (iz == z_index)
//                {
//                    topo_data_3D[iz][iy][ix] = z_m;
//                }
//                else
//                {
//                    topo_data_3D[iz][iy][ix] = z_fillValue;
//                }
//            }
        }
    }
    
    // Write data off
    
    topo_2D.putVar(&topo_data_2D[0][0]);
    
    topo_3D.putVar(&topo_data_3D[0][0][0]);
    
    delete x_data;
    delete y_data;
    delete z_data;
}

void add_national_topography(NcFile &mapfile)
{
    NcFile topography_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/oase-georef-1km-germany-2d-v01b.nc",NcFile::read);
    NcFile mapstuff_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/radolan_mapstuff.nc",NcFile::read);
    
    vector<NcDim> dimensions;
    vector<NcVar> variables;
    
    multimap<string,NcDim> topo_dims_map = topography_file.getDims();
    multimap<string,NcDim>::iterator di;
    for (di=topo_dims_map.begin(); di!=topo_dims_map.end(); ++di)
    {
        if (di->first == "t") continue;
        dimensions.push_back(di->second);
        variables.push_back(topography_file.getVar(di->first));
    }
    
    // Might use this for converting between grid
    // and projection coordinate
    CoordinateSystem<T> *cs = new CoordinateSystem<T>(dimensions, variables);
    
    // This is the main source of information
    // Note: this data is in [m].
    
    static float topography[900][900];
    for (size_t iy=0; iy < 900; iy++)
    for (size_t ix=0; ix < 900; ix++)
    topography[iy][ix]=0.0f;
    
    topography_file.getVar("topography_srtm").getVar(&topography[0][0]);
    
    // create the output file with dims/vars
    
    nc_redef(mapfile.getId());
    
    // dims
    
    NcDim dx = mapfile.addDim("national_x", 900);
    NcDim dy = mapfile.addDim("national_y", 900);
    NcDim dz = mapfile.addDim("national_z", 60);
    
    // x
    
    NcVar x = mapfile.addVar("national_x", "float", "national_x");
    x.putAtt("units", "km");
    x.putAtt("long_name", "projection_x_coordinate, national composite");
    x.putAtt("valid_min", ncFloat, -523.4622);
    x.putAtt("valid_max", ncFloat, -367.5378);
    
    // y
    
    NcVar y = mapfile.addVar("national_y", "float", "national_y");
    y.putAtt("units", "km");
    y.putAtt("long_name", "projection_y_coordinate, national composite");
    y.putAtt("valid_min", ncFloat, -4658.645f);
    y.putAtt("valid_max", ncFloat, -3759.645f);
    
    // z
    
    NcVar z = mapfile.addVar("national_z", "float", "national_z");
    z.putAtt("units", "km");
    z.putAtt("long_name", "height above sea level, national composite");
    z.putAtt("valid_min", ncFloat, 0);
    z.putAtt("valid_max", ncFloat, 14.75);
    
    double z_fillValue = -9999.0f;
    
    // Create 2D topo variable
    
    vector<NcDim> dims_2D;
    dims_2D.push_back(dy);
    dims_2D.push_back(dx);
    
    NcVar topo_2D = mapfile.addVar("national_topo_2D", ncDouble, dims_2D);
    topo_2D.putAtt("_FillValue", ncDouble, z_fillValue);
    topo_2D.putAtt("units", "m");
    topo_2D.putAtt("long_name", "height above sea level");
    topo_2D.putAtt("valid_min", ncDouble, -163.0);
    topo_2D.putAtt("valid_max", ncDouble, 14750.0f);
    topo_2D.putAtt("scale_factor", ncDouble, 1.0);
    topo_2D.putAtt("add_offset", ncDouble, 0.0);
    
    // Create 3D topo variable
    
    vector<NcDim> dims_3D;
    dims_3D.push_back(dz);
    dims_3D.push_back(dy);
    dims_3D.push_back(dx);
    
    NcVar topo_3D = mapfile.addVar("national_topo_3D", ncDouble, dims_3D);
    topo_3D.putAtt("_FillValue", ncDouble, z_fillValue);
    topo_3D.putAtt("units", "m");
    topo_3D.putAtt("long_name", "height above sea level");
    topo_3D.putAtt("valid_min", ncDouble, -163.0f);
    topo_3D.putAtt("valid_max", ncDouble, 14750.0f);
    topo_3D.putAtt("scale_factor", ncDouble, 1.0);
    topo_3D.putAtt("add_offset", ncDouble, 0.0);
    
    nc_enddef(mapfile.getId());
    
    // Topo file has no z data, therefore z must be created
    //
    float *z_data = (float *)malloc(sizeof(float)*60);
    for (size_t i=0; i<60; i++)
    {
        z_data[i] = 0.25f * ((float)i);
    }
    
    const double *y_data = cs->get_dimension_data_ptr(topography_file.getVar("y"));
    const double *x_data = cs->get_dimension_data_ptr(topography_file.getVar("x"));
    try
    {
        x.putVar(&x_data[0]);
        y.putVar(&y_data[0]);
        z.putVar(&z_data[0]);
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }
    
    // Write variables
    
    static double topo_data_2D[900][900];
    static double topo_data_3D[60][900][900];
    
    for (size_t iy=0; iy < 900; iy++)
    {
        for (size_t ix=0; ix < 900; ix++)
        {
            // Read the z value from the topography
            
            double z_m = topography[iy][ix];
            
            topo_data_2D[iy][ix] = z_m;
            
            double z_km = z_m / 1000.0f; // [m]->[km]
            
            int index_z = (int)ceil( z_km / 0.25f );
            
            if (index_z > 60)
                index_z = 60;
            else if (index_z < 0)
                index_z = 0;
            
            for (size_t iz=0; iz<60; iz++)
            {
                topo_data_3D[iz][iy][ix] = (iz<=index_z) ? z_m : z_fillValue;
            }
        }
    }
    
    // Write data off
    
    topo_2D.putVar(&topo_data_2D[0][0]);
    
    topo_3D.putVar(&topo_data_3D[0][0][0]);
    
    delete cs;
    delete z_data;
}

void add_national_mapstuff_2D(NcFile &mapfile)
{
    NcFile topography_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/oase-georef-1km-germany-2d-v01b.nc",NcFile::read);
    NcFile mapstuff_file("/Users/simon/Projects/Meteo/Ertel/data/maps/mapstuff/radolan_mapstuff.nc",NcFile::read);
    
    vector<NcDim> dimensions;
    vector<NcDim> local_dimensions;
    vector<NcVar> variables;
    
    multimap<string,NcDim> topo_dims_map = topography_file.getDims();
    multimap<string,NcDim>::iterator di;
    for (di=topo_dims_map.begin(); di!=topo_dims_map.end(); ++di)
    {
        if (di->first == "t") continue;
        dimensions.push_back(di->second);
        variables.push_back(topography_file.getVar(di->first));
    }

    local_dimensions.push_back(mapfile.getDim("national_y"));
    local_dimensions.push_back(mapfile.getDim("national_x"));

    // Might use this for converting between grid
    // and projection coordinate
    CoordinateSystem<T> *cs = new CoordinateSystem<T>(dimensions, variables);

    double z_fillValue = -9999.0f;
    
    // create the output file with dims/vars
    
    nc_redef(mapfile.getId());
    
    // Rivers
    
    NcDim rd = mapfile.addDim("riverdim", mapstuff_file.getDim("riverdim").getSize());
    NcVar rivers = mapfile.addVar("national_rivers_2D", ncFloat, local_dimensions);
    rivers.putAtt("_FillValue", ncFloat, z_fillValue);
    rivers.putAtt("units", "km");
    rivers.putAtt("long_name", "national_rivers_2D");
    rivers.putAtt("valid_min", ncFloat, 0.0);
    rivers.putAtt("valid_max", ncFloat, 1.0);
    rivers.putAtt("scale_factor", ncFloat, 1.0);
    rivers.putAtt("add_offset", ncFloat, 0.0);

    // Borders
    
    NcDim bd = mapfile.addDim("borderdim", mapstuff_file.getDim("borderdim").getSize());
    NcVar borders = mapfile.addVar("national_borders_2D", ncFloat, local_dimensions);
    borders.putAtt("_FillValue", ncFloat, z_fillValue);
    borders.putAtt("units", "km");
    borders.putAtt("long_name", "national_borders_2D");
    borders.putAtt("valid_min", ncFloat, 0.0);
    borders.putAtt("valid_max", ncFloat, 1.0);
    borders.putAtt("scale_factor", ncFloat, 1.0);
    borders.putAtt("add_offset", ncFloat, 0.0);

    nc_enddef(mapfile.getId());

    // Rivers

    float river_x[rd.getSize()];
    float river_y[rd.getSize()];
    
    mapstuff_file.getVar("river_x").getVar(&river_x[0]);
    mapstuff_file.getVar("river_y").getVar(&river_y[0]);
    
    static double river_data[900][900];
    
    for (size_t iy=0; iy < 900; iy++)
        for (size_t ix=0; ix < 900; ix++)
            river_data[iy][ix] = z_fillValue;

    typedef set<CoordinateSystem<double>::GridPoint> gp_set_t;
    
    gp_set_t line_points;
    
    for (size_t i=0; i < rd.getSize(); i++)
    {
        CoordinateSystem<double>::Coordinate coord = cs->newCoordinate();
        coord[0] = river_x[i];
        coord[1] = river_y[i];
        
        if (isnan(coord[0]) || isnan(coord[1]))
        {
            // plot the line
            
            gp_set_t::iterator li;
            
            CoordinateSystem<double>::GridPoint *lp = NULL;
            
            cout << "drawing line:" << endl;
            
            for (li=line_points.begin(); li!=line_points.end(); li++)
            {
                CoordinateSystem<double>::GridPoint p = *li;
                
                cout << p << endl;
                
                if (lp!=NULL)
                {
                    int min_x = lp->at(0) < p[0] ? lp->at(0) : p[0];
                    int max_x = lp->at(0) > p[0] ? lp->at(0) : p[0];
                    
                    int min_y = lp->at(1) < p[1] ? lp->at(1) : p[1];
                    int max_y = lp->at(1) > p[1] ? lp->at(1) : p[1];
                    
                    for (int x=min_x; x <= max_x; x++)
                    {
                        for (int y=min_y; y <= max_y; y++)
                        {
                            river_data[y][x] = 1.0;
                        }
                    }
                }

                lp = &p;
            }
            
            line_points.clear();
        }
        else
        {
            CoordinateSystem<double>::GridPoint gp = cs->newGridPoint();

            cs->reverse_lookup(coord, gp);
            
            cout << "coordinate=" << coord << " => gridpoint=" << gp << endl;
            
            line_points.insert(gp);
        }
    }
    
    // Write data off
    
    rivers.putVar(&river_data[0][0]);
    
    // Borders
    
//    float border_x[rd.getSize()];
//    float border_y[rd.getSize()];
//    
//    mapstuff_file.getVar("border_x").getVar(&border_x[0]);
//    mapstuff_file.getVar("border_y").getVar(&border_y[0]);
//    
//    static double border_data[900][900];
//    
//    for (size_t iy=0; iy < 900; iy++)
//        for (size_t ix=0; ix < 900; ix++)
//            border_data[iy][ix] = z_fillValue;
//    
//    for (size_t i=0; i < bd.getSize(); i++)
//    {
//        if (border_x[i] == std::numeric_limits<float>::quiet_NaN()
//            || border_y[i] == std::numeric_limits<float>::quiet_NaN())
//            continue;
//        
//        CoordinateSystem<double>::Coordinate coord = cs->newCoordinate();
//        coord[0] = border_x[i];
//        coord[1] = border_y[i];
//        
//        CoordinateSystem<double>::GridPoint gp = cs->newGridPoint();
//        cs->reverse_lookup(coord, gp);
//        
//        // Mark border position with 1.0
//        int ix = gp[0];
//        int iy = gp[1];
//        
//        border_data[iy][ix] = 1.0;
//    }
//    
//    // Write data off
//    
//    borders.putVar(&border_data[0][0]);

    
    delete cs;
}

void do_it()
{
    NcFile mapfile("oase-mapdata.nc",NcFile::replace,NcFile::classic);
    mapfile.putAtt("conventions","CF-1.6");
    mapfile.putAtt("authors","Jürgen Simon, Malte Diederich");
    mapfile.putAtt("created","Oct 5th 2013 18:25:00 CET");
    mapfile.putAtt("version","1.0");
    
    add_local_topography(mapfile);
    add_national_topography(mapfile);
    add_national_mapstuff_2D(mapfile);
}

/**
 *
 *
 */
int main(int argc, char** argv)
{
    using namespace cfa::utils::timer;
    using namespace cfa::utils::visit;
    using namespace m3D;
    
    // Select the correct point factory
    PointFactory<T>::set_instance( new M3DPointFactory<T>() );
    
    try
    {
        do_it();
    }
    catch (std::exception &e)
    {
        cerr << "ERROR:" << e.what() << endl;
    }
    return 0;
};
