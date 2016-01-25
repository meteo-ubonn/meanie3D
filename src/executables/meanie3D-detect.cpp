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

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/exception/info.hpp>

#include <exception>
#include <fstream>
#include <limits>
#include <locale>
#include <map>
#include <netcdf>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include <meanie3D/meanie3D.h>

using namespace std;
using namespace boost;
using namespace netCDF;
using namespace m3D;

#pragma mark -
#pragma mark Constants & Types

/** This defines the numerical type for the all variables 
 */
typedef double FS_TYPE;

#pragma mark -
#pragma mark Other command line

void print_compile_time_options() {

#if WRITE_FEATURESPACE
    cout << "\tWRITE_FEATURESPACE=1" << endl;
#else
    cout << "\tWRITE_FEATURESPACE=0" << endl;
#endif

#if WRITE_MODES
    cout << "\tWRITE_MODES=1" << endl;
#else
    cout << "\tWRITE_MODES=0" << endl;
#endif

#if WRITE_OFF_LIMITS_MASK
    cout << "\tWRITE_OFF_LIMITS_MASK=1" << endl;
#else
    cout << "\tWRITE_OFF_LIMITS_MASK=0" << endl;
#endif

#if WRITE_ZEROSHIFT_CLUSTERS
    cout << "\tWRITE_ZEROSHIFT_CLUSTERS=1" << endl;
#else
    cout << "\tWRITE_ZEROSHIFT_CLUSTERS=0" << endl;
#endif
    cout << endl;
}

#pragma mark -
#pragma mark Main

int main(int argc, char **argv) {
    using namespace m3D;

    detection_params_t<FS_TYPE> detection_params
            = Detection<FS_TYPE>::defaultParams();

    tracking_param_t tracking_params
            = Tracking<FS_TYPE>::defaultParams();

    // Declare the supported options.
    program_options::options_description desc("Options");
    utils::add_standard_options(desc);
    
    desc.add_options() ("track", "If present, tracking step is performed immediately "
            "after clustering. Required --previous-output and other "
            "clustering parameters to be set");

    add_detection_options<FS_TYPE>(desc,detection_params);
    add_tracking_options<FS_TYPE>(desc,tracking_params);
    utils::add_standard_options(desc);

    // Parse the command line
    program_options::variables_map vm;
    try {
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);
    } catch (std::exception &e) {
        cerr << "Error parsing command line: " << e.what() << endl;
        cerr << "Check meanie3D-detect --help for command line options" << endl;
        exit(EXIT_FAILURE);
    }

    // Get the command line content
    Verbosity verbosity;
    bool perform_tracking = false;
    try {
        utils::get_standard_options(vm,desc,verbosity);
        get_detection_parameters(vm,detection_params);
        utils::set_vtk_dimensions_from_args<FS_TYPE>(vm, detection_params.dimensions);
        detection_params.verbosity = verbosity;
        perform_tracking = vm.count("track") > 0;
        if (perform_tracking) {
            get_tracking_parameters<FS_TYPE>(vm,tracking_params);
            tracking_params.verbosity = verbosity;
        }
    } catch (const std::exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    // Initialize the context beforehand to allow giving the user
    // some feedback before the run.
    detection_context_t<FS_TYPE> detection_context;
    Detection<FS_TYPE>::initialiseContext(detection_params, detection_context);
    
    // Give some user feedback on the choices
    if (detection_params.verbosity > VerbositySilent) {
        cout << "----------------------------------------------------" << endl;
        cout << "meanie3D detection" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << endl;
        print_detection_params(detection_params, detection_context, vm);
        if (perform_tracking) {
            print_tracking_params(tracking_params,vm);
        }
        if (verbosity > VerbosityNormal) {
            cout << endl;
            cout << "Compiled options (use -D to switch them on and off during compile time)" << endl;
            print_compile_time_options();
        }
    }

    // Off we go
    Detection<FS_TYPE>::run(detection_params, detection_context);

    if (perform_tracking) {
        
    }

    Detection<FS_TYPE>::cleanup(detection_params, detection_context);

    return EXIT_SUCCESS;
}