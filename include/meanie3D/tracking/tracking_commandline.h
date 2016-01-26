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

#ifndef M3D_TRACKING_COMMANDLINE_H
#define	M3D_TRACKING_COMMANDLINE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/tracking/tracking.h>

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <netcdf>

namespace m3D {
    
    /**
     * Add the command line options pertaining tracking.
     * 
     * @param desc
     * @param params
     */
    template <typename T>
    void 
    add_tracking_options(program_options::options_description &desc,
                              const tracking_param_t &params) 
    {
            desc.add_options()
            ("previous", program_options::value<string>(), "Previous cluster file (netCDF)")
            ("current", program_options::value<string>(), "Current cluster file (netCDF)")
            ("tracking-variable", program_options::value<string>()->default_value(params.tracking_variable), "Variable used for histogram correlation. Must be specified when histogram weight --wt is not zero")
            ("wr", program_options::value<T>()->default_value(params.range_weight), "Weight for range correlation [0..1]")
            ("ws", program_options::value<T>()->default_value(params.size_weight), "Weight for size correlation [0..1]")
            ("wt", program_options::value<T>()->default_value(params.correlation_weight), "Weight for histogram rank correlation [0..1]")
            ("use-displacement-vectors", "If present, the algorithm uses the displacement vectors from the previous tracking result (if present) to shift clusters from the previous file to improve tracking (experimental).")
            ("merge-split-threshold", program_options::value<T>()->default_value(params.mergeSplitThreshold), "Percentage of area covered between previous/new clusters for split/merge calculation")
            ("merge-split-continuation-threshold", program_options::value<T>()->default_value(params.mergeSplitContinuationThreshold), "Minimum percentage of area covered between previous/new clusters to continue ID")
            ("discontinue-id-in-merge-and-split", "If present, the tracking discontinues cluster IDs when merging and splitting. Otherwise the largest candidate carries the ID on if the overlap exceeds --merge-split-continuation-threshold.")
            ("max-speed", program_options::value<T>()->default_value(params.maxVelocity.get()), "Maximum allowed object speed (m/s)")
            ("max-time", program_options::value<T>()->default_value(params.max_deltaT.get()), "Maximum allowed time difference between files (seconds)")
            ("max-size-deviation", 
                    program_options::value<double>()->default_value(params.max_size_deviation), 
                    "Maximum allowed difference in sizes (number of points) between clusters from previous and current file in percent [0..1]")
            #if WITH_VTK
            ("write-vtk,k", "Write out the clusters as .vtk files for visit")
            #endif
            ;
    }

    /**
     * 
     * Obtains command line parameters pertaining tracking.
     * 
     * @param vm
     * @param params
     * @param inline_tracking
     */
    template <typename T>
    void 
    get_tracking_parameters(program_options::variables_map vm, 
                            tracking_param_t &params,
                            bool inline_tracking = false)
    {
        if (vm.count("version") != 0) {
            cout << m3D::VERSION << endl;
            exit(EXIT_FAILURE);
        }

        if (! inline_tracking) {
            if (vm.count("previous") == 0) {
                cerr << "Missing previous cluster file argument --previous" << endl;
                exit(1);
            }

            try {
                params.previous_filename = vm["previous"].as<string>();
                NcFile previous_cluster_file(params.previous_filename, NcFile::read);
            } catch (const netCDF::exceptions::NcException &e) {
                cerr << "Exception opening file " << params.previous_filename << ":" << endl;
                cerr << "FATAL:" << e.what() << endl;
                exit(EXIT_FAILURE);
            }

            try {
                params.current_filename = vm["current"].as<string>();
                NcFile current_cluster_file(params.current_filename, NcFile::write);
            } catch (const netCDF::exceptions::NcException &e) {
                cerr << "FATAL:exception opening file " << params.current_filename << ":" << e.what() << endl;
                exit(EXIT_FAILURE);
            }
        }

        // max-speed
        T speed = vm["max-speed"].as<T>();
        params.maxVelocity= ::units::values::meters_per_second(speed);

        // max time
        T time = vm["max-time"].as<T>();
        params.max_deltaT= ::units::values::s(time);

        // continue ID?
        params.continueIDs = !(vm.count("discontinue-id-in-merge-and-split") > 0);
        params.useDisplacementVectors = vm.count("use-displacement-vectors") > 0;

        // merge/split continuation threshold
        params.mergeSplitContinuationThreshold = vm["merge-split-continuation-threshold"].as<T>();

        // merge/split threshold
        params.mergeSplitThreshold = vm["merge-split-threshold"].as<T>();

        // Weights
        params.range_weight = vm["wr"].as<T>();
        params.size_weight = vm["ws"].as<T>();
        params.correlation_weight = vm["wt"].as<T>();

        // tracking variable
        if (vm.count("tracking-variable") > 0) {
            params.tracking_variable = vm["tracking-variable"].as<string>();
        }

        #if WITH_VTK
        // --write-vtk
        params.write_vtk = vm.count("write-vtk") > 0;
        #endif
    }

    /**
     * Prints user feedback on tracking parameter choices to stdout.
     * 
     * @param params
     * @param vm
     */
    void 
    print_tracking_params(const tracking_param_t &params,
                          const program_options::variables_map &vm) 
    {
            cout << "\tprevious cluster file = " << params.previous_filename << endl;
            cout << "\tcurrent cluster file = " << params.current_filename << endl;
            cout << "\tvariable name for histogram comparison: " << params.tracking_variable << endl;
            cout << "\tcorrelation weights: wr=" << params.range_weight << " ws=" << params.size_weight << " wt=" << params.correlation_weight << endl;
            cout << "\tmaximum speed: " << params.maxVelocity << " [m/s]" << endl;
            cout << "\tmaximum time difference: " << params.max_deltaT << " [seconds]" << endl;
            cout << "\tmerge/split threshold: " << params.mergeSplitThreshold << endl;
            cout << "\tcontinue id in merge/split: " << (params.continueIDs ? "yes" : "no") << endl;
            cout << "\tmerge/split id continuation threshold: " << params.mergeSplitContinuationThreshold << endl;
            cout << "\tusing displacement vectors to shift previous clusters: " << (params.useDisplacementVectors ? "yes" : "no") << endl;
    #if WITH_VTK
            cout << "\twriting results out as vtk:" << (params.write_vtk ? "yes" : "no") << endl;
    #endif
            cout << endl;
    }
}

#endif	/* M3D_TRACKING_COMMANDLINE_H */

