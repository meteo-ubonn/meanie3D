#include <netcdf>
#include <iostream>

#include <meanie3D/meanie3D.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <radolan/radolan.h>

using namespace std;
using namespace boost;
using namespace Radolan;
using namespace m3D;

int main(int argc, char** argv) {
    try {
        program_options::options_description desc("Options");
        desc.add_options()
                ("help,h", "show this message")
                ("version", "print version information and exit")
                ("file,f", program_options::value<string>(), "Radolan filename or directory containing radolan scans")
                ("output-dir,o", program_options::value<string>()->default_value("."), "Path to write the results to. Defaults to current directory.")
                ("threshold,t", program_options::value<float>(), "Value threshold (depends of product)")
                ("netcdf,n", "Write scan out in netCDF/CF-Metadata format")
#if WITH_VTK
                ("vtk,k", "Write scan out in .vtk format")
#endif
                ;

        program_options::variables_map vm;
        program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
        program_options::notify(vm);

        // Version

        if (vm.count("version") != 0) {
            cout << m3D::VERSION << endl;
            exit(-1);
        }

        if (vm.count("help") == 1 || argc < 2 || vm.count("file") == 0) {
            cout << desc << "\n";
            return 1;
        }

        cout << "Endianess: " << (isLittleEndian() ? "Little" : "Big") << endl;

        // check parameters

        bool convert_to_netcdf = (vm.count("netcdf") == 0 && vm.count("vtk") == 0) || (vm.count("netcdf") > 0);

        bool convert_to_vtk = (vm.count("vtk") > 0);

        if (vm.count("file") == 0) {
            cerr << "No input" << endl;
            exit(-1);
        }

        string infile = vm["file"].as<string>();

        // File or Directory?

        std::string ending("---bin");

        namespace fs = boost::filesystem;

        fs::path path(infile);

        fs::directory_iterator end_iter;

        vector<std::string> file_paths;

        if (fs::exists(path)) {
            if (fs::is_directory(path)) {
                for (fs::directory_iterator dir_iter(path); dir_iter != end_iter; ++dir_iter) {
                    if (fs::is_regular_file(dir_iter->status())) {
                        std::string fn = dir_iter->path().generic_string();

                        if (0 == fn.compare(fn.length() - ending.length(), ending.length(), ending)) {
                            file_paths.push_back(fn);
                        }
                    }
                }
            } else {
                std::string fn = path.generic_string();
                file_paths.push_back(fn);
            }
        } else {
            cerr << "File or path does not exist: " << infile << endl;
            exit(-1);
        }

        boost::filesystem::path outpath(vm["output-dir"].as<std::string>());

        if (!boost::filesystem::exists(outpath) || !boost::filesystem::is_directory(outpath)) {
            cerr << "Can't write to path " << outpath << endl;
            exit(-1);
        }

        RDDataType *threshold = NULL;

        if (vm.count("threshold") > 0) {
            threshold = (RDDataType *) malloc(sizeof (RDDataType));
            *threshold = vm["threshold"].as<RDDataType>();
        }

        if (convert_to_netcdf) {
            vector<std::string>::iterator fi;

            for (fi = file_paths.begin(); fi != file_paths.end(); fi++) {
                std::string fn = *fi;

                boost::filesystem::path path = outpath;
                path /= boost::filesystem::path(fn).filename();
                path += ".nc";

                cout << "Converting " << fn << " to " << path.generic_string() << " ...";

                netCDF::NcFile *file = NULL;

                try {
                    file = CFConvertRadolanFile(fn.c_str(), path.generic_string().c_str(), threshold);
                    cout << " done." << endl;
                }                catch (CFFileConversionException e) {
                    cerr << endl << "Exception:" << e.what() << endl;
                }

                delete file;
            }
        }

#if WITH_VTK

        if (convert_to_vtk) {
            vector<std::string>::iterator fi;

            for (fi = file_paths.begin(); fi != file_paths.end(); fi++) {
                std::string fn = *fi;

                string outfile = fn + ".vtk";

                cout << "Converting " << fn << " to " << outfile << " ...";

                VisitUtils<double>::write_radolan_vtk(fn, outfile, threshold);

                cout << "done." << endl;
            }
        }

#endif

    }    catch (const std::exception& e) {
        cerr << "Exception: " << e.what() << endl;
    }

}
