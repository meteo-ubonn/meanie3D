#include <cf-algorithms/cf-algorithms.h>

#include <radolan/radolan.h>
#include <cf-algorithms/cf-algorithms.h>
#include <netcdfcpp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <vector>
#include <map>

typedef vector<int> coordinate_t;

/** Create gaussian distributed variable with mean 0.5 and variance
 * 
 */
float gaussian_random(int n)
{
    float sum;

    for (int i = 0; i < n; i++)
    {
        sum += rand();
    }

    return sum / ((float) n);
}

/** @returns uniform random variable (0..1)
 */
float ranf()
{
    return ((float) rand()) / ((float) RAND_MAX);
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

    if (!initialized_rand)
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
    do
    {
        x1 = 2.0 * ranf() - 1.0;
        x2 = 2.0 * ranf() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    //		use_last = 1;
    //	}

    return ( m + y1 * s);
}

/** Creates a random coordinate vector.
 * @param dimension
 * @param mean
 * @param deviation
 * @return allocated random coordinate
 */
coordinate_t *randomPoint(size_t dim, float m, float s)
{
    coordinate_t *newPoint = new coordinate_t(dim);

    for (size_t i = 0; i < dim; i++)
    {
        newPoint->at(i) = box_muller(m, s);
    }

    return newPoint;
}

/** Write out an axis with dim->num_dim() grid points, values varying
 * linear between min and max
 */
void writeAxis(NcFile* file, NcDim *dim, float min, float max)
{
    // write x-axis information. Make 0 the center of the map. 

    float *data = (float *) malloc(sizeof (float) * dim->size());

    for (int i = 0; i < dim->size(); i++)
    {
        data[i] = min + ((max - min) / dim->size()) * i;
    }

    NcVar *dimData = file->add_var(dim->name(), ncFloat, dim);

    dimData->add_att("value_min", min);

    dimData->add_att("value_max", max);

    dimData->put(data, dim->size());

    free(data);
}

/**
 */
void writeCloud(NcFile* file, NcVar* variable, vector<NcDim *> dims, vector<float> center, vector<float> mean, vector<float> deviation, size_t cloudSize)
{
    map< NcDim*, NcVar* > vars;

    for (size_t i = 0; i < dims.size(); i++)
    {
        vars[ dims.at(i) ] = file->get_var(dims.at(i)->name());
    }

    // writing individual points here, counts are all 1

    long *counts = (long *) malloc(sizeof (long) * dims.size());

    for (size_t i = 0; i < dims.size(); i++)
    {
        counts[i] = 1;
    }

    // allocate a cursor

    long *cursor = (long *) malloc(sizeof (long) * dims.size());

    // start generating random points (with random values between 0 and 1)

    size_t numPoints = 0;

    do
    {
        // genrate a random coordinate

        for (size_t d = 0; d < dims.size(); d++)
        {
            bool valid = false;

            while (!valid)
            {
                NcDim *dim = dims.at(d);

                float rand = box_muller(mean.at(d), deviation.at(d));

                // re-transform to grid coordinates

                float min = vars[dim]->get_att("value_min")->values()->as_double(0);

                float max = vars[dim]->get_att("value_max")->values()->as_double(0);

                long n = (long) round((dim->size() - 1)*(rand - min) / (max - min));

                if (n >= 0 && n < dim->size())
                {
                    cursor[d] = n;

                    valid = true;
                }
            }
        }

        // generate a random value

        float value[] = {1.0f}; //ranf();

        variable->set_cur(cursor);

        variable->put(value, counts);

        numPoints++;

    } while (numPoints < cloudSize);

    free(counts);

    free(cursor);
}

void writeCloudND(const char *filename, size_t cloud_size, size_t gridSize, float mean, float deviation, vector<float> *centerOffset = NULL) {
 }

void writeCloud2D(const char *filename, size_t cloud_size, size_t gridSize, float m, float s, vector<float> *centerOffset = NULL)
{
    NcFile *file = new NcFile(filename, NcFile::Replace);

    // Create gaussian variable with 3 dimensions
    vector<NcDim *> dims;

    NcDim* x = file->add_dim("x", gridSize);
    writeAxis(file, x, -100.0, 100.0);
    dims.push_back(x);

    NcDim* y = file->add_dim("y", gridSize);
    writeAxis(file, y, -100.0, 100.0);
    dims.push_back(y);

    NcVar *gaussian = file->add_var("gaussian", ncFloat, x, y);
    gaussian->add_att("valid_min", 0.0);
    gaussian->add_att("valid_max", 1.0);

    // write cloud out

    vector<float> mean;
    mean.push_back(m);
    mean.push_back(m);

    vector<float> deviation;
    deviation.push_back(s);
    deviation.push_back(s);

    vector<float> center;

    if (centerOffset != NULL)
    {
        center = *centerOffset;
    } else
    {
        center.push_back(0.0);
        center.push_back(0.0);
    }

    writeCloud(file, gaussian, dims, center, mean, deviation, cloud_size);

    file->close();
}

void writeCloud3D(const char *filename, size_t cloud_size, size_t gridSize, float m, float s, vector<float> *centerOffset = NULL)
{
    NcFile *file = new NcFile(filename, NcFile::Replace);

    // Create gaussian variable with 3 dimensions
    vector<NcDim *> dims;

    NcDim* x = file->add_dim("x", gridSize);
    writeAxis(file, x, -100.0, 100.0);
    dims.push_back(x);

    NcDim* y = file->add_dim("y", gridSize);
    writeAxis(file, y, -100.0, 100.0);
    dims.push_back(y);

    NcDim* z = file->add_dim("z", gridSize);
    writeAxis(file, z, -100.0, 100.0);
    dims.push_back(z);

    NcVar *gaussian = file->add_var("gaussian", ncFloat, x, y, z);
    gaussian->add_att("valid_min", 0.0);
    gaussian->add_att("valid_max", 1.0);

    // write cloud out

    vector<float> mean;
    mean.push_back(m);
    mean.push_back(m);
    mean.push_back(m);

    vector<float> deviation;
    deviation.push_back(s);
    deviation.push_back(s);
    deviation.push_back(s);

    vector<float> center;

    if (centerOffset != NULL)
    {
        center = *centerOffset;
    } else
    {
        center.push_back(0.0);
        center.push_back(0.0);
        center.push_back(0.0);
    }

    writeCloud(file, gaussian, dims, center, mean, deviation, cloud_size);

    file->close();
}

int main(int argc, char** argv)
{
    writeCloud2D("gaussian2D.nc", 2000, 101, 0.0, 40.0);

    // writeCloud3D("gaussian3D.nc", 4000, 101, 0.0, 40.0 );

    return 0;
}