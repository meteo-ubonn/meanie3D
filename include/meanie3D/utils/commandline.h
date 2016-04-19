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

#ifndef M3D_COMMANDLINE_H
#define	M3D_COMMANDLINE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils/visit.h>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include <string>
#include <vector>

namespace m3D {
    namespace utils {

        using namespace boost;

        /**
         * Parses a command line parameter with comma-separated values
         * into a vector of strings.
         * 
         * @param vm
         * @param name
         * @return 
         */
        std::vector<std::string>
        parse_string_vector(const program_options::variables_map &vm,
                const std::string &name) 
        {
            if (vm.count(name) == 0) {
                cerr << "FATAL: could not find parameter '" << name << "'" << endl;
                exit(EXIT_FAILURE);
            }
            
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            std::vector<std::string> result;
            std::string value = vm[name].as<std::string>();
            tokenizer tokens(value, char_separator<char>(","));
            tokenizer::iterator ti;
            for (ti = tokens.begin(); ti != tokens.end(); ++ti) {
                result.push_back(*ti);
            }
            return result;
        }
        
        /**
         * 
         * @param desc
         */
        void
        add_standard_options(program_options::options_description &desc, 
                bool with_vtk_dimensions = true)  {
            desc.add_options()
                ("help,h", "produce help message")
                ("version", "print version information and exit")
                ("verbosity", 
                    program_options::value<unsigned int>()->default_value(1u), 
                    "Verbosity level [0..3]. (0=silent, 1=normal, 2=show details, 3=show all details)");
            #if WITH_VTK
                if (with_vtk_dimensions) {
                    desc.add_options()
                        ("vtk-dimensions", 
                            program_options::value<string>(),
                            "VTK files are written in the order of dimensions given. This may lead to wrong results if the order of the dimensions is not x,y,z. Add the comma-separated list of dimensions here, in the order you would like them to be written as (x,y,z)");
                }
            #endif
        }
        
        /**
         * Handles some common parameters, like:
         * <ul>
         *  <li>--help</li>
         *  <li>--version</li>
         *  <li>--verbosity</li>
         * </ul>
         * @param vm
         * @param desc
         * @param verbosity
         */
        void
        get_standard_options(const int &argc,
                             const program_options::variables_map &vm,
                             const program_options::options_description &desc,
                             Verbosity &verbosity) 
        {
            // --help
            if (vm.count("help") == 1 || argc < 2) {
                cout << desc << "\n";
                exit(EXIT_FAILURE);
            }

            // --version
            if (vm.count("version") == 1) {
                cout << ::m3D::VERSION << endl;
                exit(EXIT_SUCCESS);
            }

            // --verbosity
            unsigned int vb = vm["verbosity"].as <unsigned int> ();
            if (vb > VerbosityAll) {
                cerr << "Illegal value for parameter --verbosity. Only values from 0 .. 3 are allowed" << endl;
                exit(EXIT_FAILURE);
            } else {
                verbosity = (Verbosity) vb;
            }
        }

       
        /**
         * Checks if --vtk-dimensions is present. If yes, creates and
         * sets a dimension mapping for visit utils.
         * 
         * @param vm command line
         * @param dimensions dimension names to map to
         */
        template <typename T>
        void
        set_vtk_dimensions_from_args(const program_options::variables_map &vm,
                const std::vector<std::string> &dimensions) 
        {
            #if WITH_VTK
            // --vtk-dimensions
            if (vm.count("vtk-dimensions") > 0) {
                vector<std::string> vtk_dimensions = parse_string_vector(vm,"vtk-dimensions");
                vector<size_t> vtk_dimension_indexes;
                for (size_t i=0; i < vtk_dimensions.size(); i++) {
                    bool found = false;
                    for (size_t j=0; j < dimensions.size() && !found; j++) {
                        if (vtk_dimensions[i] == dimensions[j]) {
                            found = true;
                            vtk_dimension_indexes.push_back(j);
                            found = true;
                        }
                    }
                    if (!found) {
                        cerr << "FATAL:could not map dimensions. "
                             << "Please check --vtk-dimensions" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                VisitUtils<T>::VTK_DIMENSION_INDEXES = vtk_dimension_indexes;
            }
            #endif
        }
    }
}

#endif	/* M3D_COMMANDLINE_H */

