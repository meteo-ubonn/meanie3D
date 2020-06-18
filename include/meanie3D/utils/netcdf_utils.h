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

#ifndef M3D_NETCDFUTILS_H
#define M3D_NETCDFUTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/cast.hpp>

#include <vector>
#include <netcdf>
#include <iostream>
#include <exception>
#include <string.h>

namespace m3D {
    namespace utils {
        namespace netcdf {

            using namespace netCDF;

            /** Figures out the number of values in the given variable by multiplying
             * all it's dimension's sizes.
             * @param variable
             * @return number of values
             */
            size_t num_vals(const NcVar &var) {
                size_t count = 1;
                for (int dim_index = 0; dim_index < var.getDimCount(); dim_index++) {
                    count = count * var.getDim(dim_index).getSize();
                }
                return count;
            }

            template<typename T>
            T *
            readNcByte(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                unsigned short *values = (unsigned short *) calloc(num_vals(variable), sizeof(unsigned short));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            template<typename T>
            T *readNcChar(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                char *values = (char *) calloc(num_vals(variable), sizeof(char));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            template<typename T>
            T *readNcDouble(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                double *values = (double *) calloc(num_vals(variable), sizeof(double));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            template<typename T>
            T *readNcFloat(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                float *values = (float *) calloc(num_vals(variable), sizeof(float));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            template<typename T>
            T *readNcInt(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                long *values = (long *) calloc(num_vals(variable), sizeof(long));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            template<typename T>
            T *readNcShort(const NcVar &variable, const T &scale_factor, const T &add_offset) {
                short *values = (short *) calloc(num_vals(variable), sizeof(short));
                variable.getVar(values);
                T *returnValue = (T *) calloc(num_vals(variable), sizeof(T));
                for (int index = 0; index < num_vals(variable); index++) {
                    returnValue[index] = (T) values[index] * scale_factor + add_offset;
                }
                free(values);
                return returnValue;
            }

            /** Read the whole variable content into one sequential array
             * @param variable
             * @param allocated memory buffer with content
             */
            template<typename T>
            T *
            readNetCDFVariable(const NcVar &variable) {
                using namespace netCDF;
                using namespace std;

                T *returnValue = NULL;
                T scale_factor = 1.0;
                try {
                    map<std::string, NcVarAtt> attributes = variable.getAtts();
                    map<std::string, NcVarAtt>::iterator fi;
                    fi = attributes.find("scale_factor");
                    if (fi != attributes.end()) {
                        fi->second.getValues(&scale_factor);
                    }
                    T offset = 0.0;
                    fi = attributes.find("add_offset");
                    if (fi != attributes.end()) {
                        fi->second.getValues(&offset);
                    }
                    switch (variable.getType().getTypeClass()) {
                        case netCDF::NcType::nc_BYTE:
                            returnValue = readNcByte<T>(variable, scale_factor, offset);
                            break;

                        case netCDF::NcType::nc_CHAR:
                            returnValue = readNcChar<T>(variable, scale_factor, offset);
                            break;

                        case netCDF::NcType::nc_DOUBLE:
                            returnValue = readNcDouble<T>(variable, scale_factor, offset);
                            break;

                        case netCDF::NcType::nc_FLOAT:
                            returnValue = readNcFloat<T>(variable, scale_factor, offset);
                            break;

                        case netCDF::NcType::nc_INT:
                            returnValue = readNcInt<T>(variable, scale_factor, offset);
                            break;

                        case netCDF::NcType::nc_SHORT:
                            returnValue = readNcShort<T>(variable, scale_factor, offset);

                        default: {
                            std::cerr << "ERROR:variable " << variable.getName() << " has unsupported type "
                                      << variable.getType().getTypeClassName() << std::endl;
                            returnValue = NULL;
                        }
                    }
                } catch (netCDF::exceptions::NcException &e) {
                    cerr << "ERROR: can not get attributes of variable " << variable.getName() << e.what() << endl;
                }
                return returnValue;
            }

            /** Retrieves a valus from a specific gridpoint of a netcdf variable.
             * As stated by http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html#packed-data
             * the valid_min/valid_max are applied against the PACKED data, that is before scaling and offset.
             * @param variable
             * @param grid point
             * @param time_index if -1 time is ignored. If 0 or greater, 
             * the actual point is read from that index in time
             * @param reference to bool variable, which will indicate if the 
             * value is within valid range or not after the call
             * @param scale factor (for unpacking, defaults to 1)
             * @param offset (for unpacking, defaults to 0)
             * @param valid_min (defaults to the min value for the template parameter)
             * @param valid_max (defaults to the max value for the template parameter)
             */
            template<typename T>
            inline
            T retrieveValueAt(const NcVar &var,
                              const vector<int> &gridpoint,
                              int time_index,
                              bool &is_valid,
                              const T scale_factor = 1.0,
                              const T offset = 0.0,
                              const T *fill_value = NULL,
                              const T valid_min = std::numeric_limits<T>::min(),
                              const T valid_max = std::numeric_limits<T>::max()) {
                T value = 0.0;
                vector<size_t> index(gridpoint.begin(), gridpoint.end());
                if (time_index >= 0) {
                    index.insert(index.begin(), boost::numeric_cast<size_t>(time_index));
                }
                switch (var.getType().getTypeClass()) {
                    case netCDF::NcType::nc_BYTE: {
                        unsigned char val;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    case netCDF::NcType::nc_CHAR: {
                        char val;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    case netCDF::NcType::nc_DOUBLE: {
                        double val = 0.0;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    case netCDF::NcType::nc_FLOAT: {
                        float val = 0.0;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    case netCDF::NcType::nc_INT: {
                        int val = 0;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    case netCDF::NcType::nc_SHORT: {
                        short val = 0;
                        var.getVar(index, &val);
                        value = boost::numeric_cast<T>(val);
                    }
                        break;

                    default: {
                        var.getVar(index, &value);
                    }
                }

                if (fill_value != NULL) {
                    is_valid = (value != *fill_value) && (value >= valid_min && value <= valid_max);
                } else {
                    // HOTFIX for Malte's packed/unpacked issue with
                    // type byte. Remember to remove this when the data
                    // has been fixed!
                    if (var.getType().getTypeClass() == netCDF::NcType::nc_BYTE) {
                        is_valid = (value >= valid_min && value <= (valid_max / scale_factor));
                    } else {
                        is_valid = (value >= valid_min && value <= valid_max);
                    }
                }

                // scale first, then offset
                T unpacked_value = scale_factor * value + offset;
                return unpacked_value;
            }

            /**
             * Retrieves the valid_min/valid_max. Looks at valid_range,
             * valid_min and valid_max. Values are <b>packed</b>.
             *
             * @param variable
             * @param valid_min
             * @param valid_min
             * @throw std::runtime_error if neither valid_min+valid_max nor valid_range existed
             */
            template<typename T>
            void
            get_valid_range(const NcVar &var, T &valid_min, T &valid_max) {
                bool have_valid_range = false;

                // check for 'valid_range' attribute first
                T values[2];
                try {
                    var.getAtt("valid_range").getValues(&values[0]);
                    valid_min = values[0];
                    valid_max = values[1];
                    have_valid_range = true;
                } catch (const std::exception &e) {
                    have_valid_range = false;
                }

                if (!have_valid_range) {
                    // need both valid_min and valid_max
                    try {
                        var.getAtt("valid_min").getValues(&valid_min);
                        have_valid_range = true;
                    } catch (const std::exception &e) {
                        have_valid_range = false;
                    }
                    try {
                        var.getAtt("valid_max").getValues(&valid_max);
                        have_valid_range = true;
                    } catch (const std::exception &e) {
                        have_valid_range = false;
                    }
                }

                // No valid_range or valid_min+valid_max constitutes
                // a reason for exception
                if (!have_valid_range) {
                    std::stringstream msgBuf;
                    msgBuf << "Variable " << var.getName() << " does not have valid range";
                    throw std::runtime_error(msgBuf.str());
                }
            }

            /**
             * Retrieves the valid_min/valid_max. Looks at valid_range,
             * valid_min and valid_max. Values are <b>unpacked</b>.
             *
             * @param valid_min
             * @param valid_max
             * @throw std::exception if neither valid_min+valid_max nor valid_range existed
             */
            template<typename T>
            void
            unpacked_limits(const NcVar &var, T &valid_min, T &valid_max) {
                NcVarAtt att;

                T scale_factor = 1.0;
                try {
                    att = var.getAtt("scale_factor");
                    att.getValues(&scale_factor);
                } catch (netCDF::exceptions::NcAttExists &e) {
                } catch (std::exception &e) {
                }

                T offset = 0.0;
                try {
                    att = var.getAtt("add_offset");
                    att.getValues(&offset);
                } catch (netCDF::exceptions::NcAttExists &e) {
                } catch (std::exception &e) {
                }

                ::m3D::utils::netcdf::get_valid_range(var, valid_min, valid_max);
                valid_min = offset + scale_factor * valid_min;
                valid_max = offset + scale_factor * valid_max;
            }

            /** Retrieves an simple attribute value with some checks
             * @param netcdf variable
             * @param attribute name
             * @return value
             */
            template<typename T>
            T get_attribute_value(NcVar variable, const char *attributeName) {
                T result;

                if (variable.isNull()) {
                    std::string msg = "Can't access variable '" + variable.getName() + "'";
                    throw new invalid_argument(msg);
                }

                try {
                    NcVarAtt att = variable.getAtt(attributeName);
                    if (att.isNull()) {
                        std::string msg = "Can't access attribute '" + std::string(attributeName) + "' in variable '" +
                                          variable.getName() + "'";
                        throw new invalid_argument(msg);
                    }

                    att.getValues(&result);
                } catch (netCDF::exceptions::NcException &e) {
                    std::string msg = "Can't access attribute '" + std::string(attributeName) + "' in variable '" +
                                      variable.getName() + "'";
                    throw new invalid_argument(msg);
                }

                return result;
            }

            /**
             * Evaluates the presence of time(time). It does this by
             * evaluating the time variable's attributes. If the attribute
             * 'long_name' is present and it's value is 'time', it is 
             * assumed that the time(time) is present (according to the 
             * cf-metadata convention). If not, the weaker criterion is
             * used that the variable has one dimension and that this 
             * dimension is called either 't' or 'time'.
             *
             * @param netcdf filename
             * @param name of time variable (default 'time')
             * @return <code>true</code> if file has time(time). 
             */
            bool
            have_time(const std::string &filename,
                      const std::string &time_var_name = "time") {
                bool result = false;
                NcFile *file = NULL;
                try {
                    file = new NcFile(filename, NcFile::read);
                    if (!file->isNull()) {
                        NcVar timeVar = file->getVar(time_var_name);
                        if (!timeVar.isNull()) {

                            // check if the time dimension is explicity denoted
                            // 'time' by long name
                            try {
                                std::string long_name;
                                NcVarAtt ncLongName = timeVar.getAtt("long_name");
                                ncLongName.getValues(long_name);
                                result = (long_name == "time");
                            } catch (netCDF::exceptions::NcException &e) {
                            }

                            // If that was not the case, we can check if
                            // the time dimension is called t or time
                            if (!result) {
                                if (timeVar.getDimCount() == 1) {
                                    NcDim t = timeVar.getDim(0);
                                    result = t.getName() == "t" || t.getName() == "time";
                                }
                            }
                        }
                    }
                } catch (netCDF::exceptions::NcException &e) {
                }
                if (file != NULL) delete file;
                return result;
            }

            /** 
             * Adds a dimension and variable time eg. time(time).
             * Time is a 1-D array with the given value
             */
            bool
            add_time(NcFile *file, unsigned long timestamp,
                     bool toggle_defmode = true) {
                try {
                    if (toggle_defmode)
                        nc_redef(file->getId());
                    // Add a dimension 'time'
                    NcDim dTime = file->addDim("time", 1);
                    // Add variable 'time(time)'
                    // the NetCDF API is HORRIFIC!! Can't use anything
                    // but "int" at this point. Puking out my breakfast.
                    NcVar vTime = file->addVar("time", "int", "time");
                    vTime.putAtt("long_name", "time");
                    vTime.putAtt("units", "seconds since Jan 1 1970 00:00:00 GMT");
                    if (toggle_defmode)
                        nc_enddef(file->getId());

                    // Add value
                    unsigned long values[1] = {timestamp};
                    vTime.putVar(&values[0]);
                } catch (const netCDF::exceptions::NcException &e) {
                    cerr << "ERROR:could not add variable 'time(time)' to file '" << file->getName() << "' : "
                         << e.what() << endl;
                    return false;
                } catch (const std::exception &e) {
                    cerr << "ERROR:could not add variable 'time(time)' to file '" << file->getName() << "' : "
                         << e.what() << endl;
                    return false;
                }

                return true;
            }

            /** Adds a dimension and variable time eg. time(time).
             * Time is a 1-D array with the given value
             */
            bool add_time(std::string fn, unsigned long timestamp, bool toggle_defmode = true) {
                NcFile *file = new NcFile(fn, NcFile::write);
                bool result = add_time(file, timestamp, toggle_defmode);
                delete file;
                return result;
            }

            /** Tries to figure out time dimension and time variable. 
             * The time dimension is tried as 'time' first, then 't'. 
             * The variable is the first variable with standard_name = "time" 
             * and of sole dimension being the time dimension.
             *
             * @throws std::runtime_exception if time can not be defined
             */
            void get_time_dim_and_var(NcFile &file,
                                      NcDim &time_dim, NcVar &time_var)
            throw(std::runtime_error) {
                // find time dimension 'time' or 't'
                time_dim = file.getDim("time");
                if (time_dim.isNull()) {
                    time_dim = file.getDim("t");
                }
                if (time_dim.isNull()) {
                    throw runtime_error("ERROR:can't read time dimension");
                }

                // Find variable time(time)
                std::multimap<std::string, NcVar>::iterator vi;
                std::multimap<std::string, NcVar> vars = file.getVars();
                for (vi = vars.begin(); vi != vars.end(); ++vi) {
                    try {
                        // check if the variable depends on dimension time_dim alone
                        vector<NcDim> dims = vi->second.getDims();
                        if (dims.size() != 1) continue;
                        if (dims.at(0).getId() == time_dim.getId()) {
                            // Check for standard_name = "time"
                            try {
                                NcVarAtt standard_name = vi->second.getAtt("standard_name");
                                std::string value;
                                standard_name.getValues(value);
                                if (strcmp(value.c_str(), "time") == 0) {
                                    time_var = vi->second;
                                    break;
                                }
                            } catch (::netCDF::exceptions::NcException &e) {
                            }

                            // Check for long_name = "time"
                            try {
                                NcVarAtt long_name = vi->second.getAtt("long_name");
                                std::string value;
                                long_name.getValues(value);
                                if (strcmp(value.c_str(), "time") == 0) {
                                    time_var = vi->second;
                                    break;
                                }
                            } catch (::netCDF::exceptions::NcException &e) {
                            }

                        }
                    } catch (const ::netCDF::exceptions::NcException &e) {
                    }
                }
                if (time_var.isNull())
                    throw runtime_error("ERROR:can't read time variable");
            }

            /** 
             * Get the time value for the given index. This assumes
             * that there is a time(time) constellation in the file.
             *
             * TODO: time unit handling
             *
             * @param file
             * @param index in time dimension
             * @param name of time variable (defaults to 'time')
             *
             * @return actual value for time at time_index
             */
            template<typename T>
            T
            get_time(const std::string &filename,
                     const int &time_index,
                     const std::string &time_var_name = "time")
            throw(std::runtime_error) {
                try {
                    NcFile file(filename, NcFile::read);
                    if (file.isNull()) {
                        throw runtime_error("ERROR:can't open file " + filename);
                    }
                    NcDim time_dim;
                    NcVar time_var;
                    try {
                        get_time_dim_and_var(file, time_dim, time_var);
                    } catch (std::exception &) {
                    }

                    if (time_var.isNull() || time_dim.isNull()) {
                        // No time(time)
                        return 0;
                    } else {
                        if (time_index < time_dim.getSize() && time_index >= 0) {
                            // Have time(time) and time_index
                            T *times = new T[time_dim.getSize()];
                            time_var.getVar(times);
                            T result = times[time_index];
                            delete[] times;
                            return result;
                        } else if (time_index == NO_TIME) {
                            // Have time(time) but no time_index 
                            T time;
                            time_var.getVar(&time);
                            return time;
                        } else {
                            // Might have time(time) but index is out of range
                            return 0;
                        }
                    }
                } catch (const ::netCDF::exceptions::NcException &e) {
                    throw runtime_error("ERROR:can't access file " + filename + ":" + e.what());
                }
            }


            /**
             * Performs a checked method of retrieving time information
             * from a given file. This is a little complicated, as there
             * can be a number of cases:
             * 
             * a) time_index = NO_TIME and no time(time) exists: really no time
             * b) time_index = NO_TIME and time(time) exists: t = time(0)
             * c) time_index > NO_TIME t = time(time_index)
             * d) filename is "rico.out.xy.*" => use special treatment
             * 
             * @param filename
             * @param time_index 
             * @param time_var_name (defaults to 'time')
             * @return 
             */
            template<typename T>
            T
            get_time_checked(const std::string &filename,
                             const int &time_index,
                             const std::string &time_var_name = "time") {
                T timestamp = -1;

                bool have_t = have_time(filename, time_var_name);

                if (boost::contains(filename, "rico.out.xy.")) {

                    // TODO: this block was put in for a specific tracking
                    // inter-comparison problem. Remove when the project is through!!
                    // time for this is days since simulation start

                    // Days in simulation -> seconds
                    double time_in_days = get_time<double>(filename, time_index);
                    // a day has 24 * 60 * 60 seconds
                    timestamp = (T) round(time_in_days * 24.0 * 60.0 * 60.0);
                } else if (time_index == NO_TIME && !have_t) {
                    timestamp = 0l;
                } else if (time_index == NO_TIME && have_t) {
                    timestamp = get_time<T>(filename, 0);
                } else if (time_index != NO_TIME && have_t) {
                    timestamp = get_time<T>(filename, time_index);
                }

                return timestamp;
            }


            /** Extracts the names of the given variables and returns
             * them as vectors.
             * @param vector of netcdf variables
             * @return vector of strings
             */
            std::vector<std::string>
            to_names(const std::vector<netCDF::NcVar> &variables) {
                std::vector<std::string> result;
                for (size_t i = 0; i < variables.size(); i++) {
                    result.push_back(variables[i].getName());
                }
                return result;
            }

            /** Extracts the names of the given variables and returns
             * them as vectors.
             * @param NetCDF file
             * @param vector of strings
             * @return vector of netcdf variables
             */
            std::vector<netCDF::NcVar>
            to_variables(const NcFile *file, const std::vector<std::string> &names) {
                std::vector<netCDF::NcVar> variables;
                for (size_t i = 0; i < names.size(); i++) {
                    variables.push_back(file->getVar(names[i]));
                }
                return variables;
            }

            /** Create an NcFile copying dimensions and dimension variables
             * from the given file. Note: at this time only supports simple
             * dimension variables, aka x(x) etc.
             * 
             * @param fn
             * @return 
             */
            template<typename T>
            NcFile *
            create_file_by_copying_dimensions(const std::string &source_path,
                                              const std::string &dest_path,
                                              const std::vector<std::string> &dimensionNames) {
                NcFile *source = NULL;

                // open the original

                try {
                    source = new NcFile(source_path.c_str(), NcFile::read);
                } catch (const netCDF::exceptions::NcException &e) {
                    cerr << "FATAL:exception opening file " << source_path << " for reading : " << e.what() << endl;
                    exit(EXIT_FAILURE);
                }

                // get source file dimensions, using coordinate system
                CoordinateSystem<T> *source_cs = new CoordinateSystem<T>(source, dimensionNames);

                // open the copy

                NcFile *dest = NULL;

                try {
                    dest = new NcFile(dest_path, NcFile::write);
                } catch (const netCDF::exceptions::NcException &e) {
                    cerr << "FATAL:Exception opening file " << dest_path << " for writing : " << e.what() << endl;
                    exit(EXIT_FAILURE);
                }

                // copy dimensions and dimension data

                for (size_t di = 0; di < source_cs->rank(); di++) {
                    // dimension

                    NcDim d = source_cs->dimensions()[di];
                    NcDim dest_dim = dest->addDim(d.getName(), d.getSize());

                    // dimension variable

                    NcVar v = source_cs->dimension_variables()[di];
                    NcVar source_var = dest->addVar(v.getName(), v.getType(), dest_dim);

                    // copy dimension data

                    T *data = source_cs->get_dimension_data_ptr(di);
                    source_var.putVar(data);
                }

                // close source file

                delete source;
                delete source_cs;

                return dest;
            }

            /**
             * Copy the given list of netCDF dimensions from source file
             * to destination file.
             * 
             * @param dimensions
             * @param source
             * @param dest
             * @return 
             */
            vector <NcDim>
            copy_dimensions(const vector <string> &dimensions,
                            const NcFile *source, NcFile *dest) {
                vector<NcDim> ncDimensions;
                for (size_t di = 0; di < dimensions.size(); di++) {
                    string dimension = dimensions[di];
                    NcDim ncSourceDim = source->getDim(dimension);
                    int size = ncSourceDim.getSize();
                    NcDim ncDimension = dest->addDim(dimension, size);
                    ncDimensions.push_back(ncDimension);
                }
                return ncDimensions;
            }

            /**
             * Creates a copy of a netCDF variable in another file. This
             * operation requires, that the dimensions for the variable
             * to be copied are already in place. 
             * 
             * @param variable
             * @param source
             * @param dest
             * @param with_data Only if <code>true</code> the actual data is copied.
             * @return <code>true</code> if operation was success.
             * @throws netCDF::exceptions::NcException
             */
            template<typename T>
            bool copy_variable(const std::string &variable,
                               const NcFile *source, NcFile *dest, bool with_data) {
                // This should no longer be required with netCDF C++ 4 API
                // nc_redef(file->getId());                

                NcVar sourceVar = source->getVar(variable);

                // List of dimensions for the copy (in dest file)
                vector<NcDim> copyDims;

                // collect the dimension info
                vector<NcDim> sourceDims = sourceVar.getDims();
                for (size_t j = 0; j < sourceDims.size(); j++) {
                    NcDim sourceDim = sourceDims[j];
                    NcDim destDim = dest->getDim(sourceDim.getName());
                    if (destDim.isNull()) {
                        cerr << "ERROR:could not find dimension '"
                             << sourceDim.getName() << "'" << endl;
                        return false;
                    }
                    copyDims.push_back(destDim);
                }

                // create a dummy variable
                NcVar copy;
                copy = dest->addVar(variable, sourceVar.getType(), copyDims);
                if (copy.isNull()) {
                    cerr << "ERROR:could not create variable '"
                         << variable << "'" << endl;
                    return false;
                }
                copy.setCompression(false, true, 3);

                // Copy attributes
                map<string, NcVarAtt> attributes = sourceVar.getAtts();
                map<string, NcVarAtt>::iterator at;
                for (at = attributes.begin(); at != attributes.end(); at++) {
                    NcVarAtt a = at->second;
                    size_t size = a.getAttLength();
                    void *data = (void *) malloc(size);
                    a.getValues(data);
                    NcVarAtt att_copy = copy.putAtt(a.getName(), a.getType(), size, data);
                    free(data);

                    if (att_copy.isNull()) {
                        cerr << "ERROR: could not write attribute " << a.getName() << endl;
                        return false;
                    }
                }

                // Copy data is necessary
                if (with_data) {
                    size_t number_of_values = num_vals(sourceVar);
                    T *data = (T *) malloc(sizeof(T) * number_of_values);
                    if (data == NULL) {
                        cerr << "ERROR:out of memory" << endl;
                        return false;
                    }
                    sourceVar.getVar(data);
                    copy.putVar(data);
                    free(data);
                }

                return true;
            }
        }
    }
}

#endif
