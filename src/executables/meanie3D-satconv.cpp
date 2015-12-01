//
//  meanie3D-satconv
//

#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

using namespace std;
using namespace boost;

#pragma mark -
#pragma mark Definitions

typedef enum {
    RadicanceToTemperature,
    TemperatureToRadiance
} ConversionType;

// Radiation constant #1
const double c1 = 1.19104e-05;
// Radiation constant #2
const double c2 = 1.43877;

// Names of the variables available for conversion
const string variables[8] = {
    "msevi_l15_ir_039", 
    "msevi_l15_ir_087", 
    "msevi_l15_ir_097", 
    "msevi_l15_ir_108", 
    "msevi_l15_ir_120", 
    "msevi_l15_ir_134",
    "msevi_l15_wv_062",
    "msevi_l15_wv_073"
};

// Calibration constant alpha
const double alpha[8] = {
    0.9954, 
    0.9996, 
    0.9999,
    0.9983, 
    0.9988, 
    0.9981,
    0.9963,
    0.9991
};

// Calibration constant beta
const double beta[8] = {
    3.438, 
    0.179, 
    0.056, 
    0.64, 
    0.408, 
    0.561,
    2.185,
    0.47
};

// Wavenumber 
const double wavenum[8] = {
    2568.832, 
    1148.62, 
    1035.289, 
    931.7, 
    836.445, 
    751.792,
    1600.548,
    1360.33
};
       
#pragma mark -
#pragma mark Command line parsing

void parse_commmandline(program_options::variables_map vm, 
        int& variableIndex,
        double &value,
        ConversionType &conversionType) 
{
    if (vm.count("temperature") != 0) {
        conversionType = RadicanceToTemperature;
    } else if (vm.count("radiance") != 0) {
        conversionType = TemperatureToRadiance;
    } else {
        cerr << "Specify -r or -t" << endl;
        exit(EXIT_FAILURE);
    }
    
    if (vm.count("value") == 0) {
        cerr << "Missing --value argument" << endl;
        exit(EXIT_FAILURE);
    }
    value = vm["value"].as<double>();
    
    if (vm.count("variable") == 0) {
        cerr << "Missing --variable argument" << endl;
        exit(EXIT_FAILURE);
    }
    string variable = vm["variable"].as<string>();
    
    variableIndex = -1;
    for (int i=0; i<8; i++) {
        if (variables[i] == variable) {
            variableIndex = i;
            break;
        }
    }
    if (variableIndex < 0) {
        cerr << "Illegal value for --variable" << endl;
        exit(EXIT_FAILURE);
    }

}

#pragma mark -
#pragma mark Helper Methods
        
/** Calculates the brightness temperature in degree centigrade
* from the given spectral radiance.
 * 
* @param index of the variable to be converted
* @param radiance value for the given channel
* @return equivalent brightness temperature in [degree C]
*/
double brightness_temperature(const int var_index, const double &radiance) {
   double nu = wavenum[var_index];
   double Tbb = c2 * nu / log(1 + nu * nu * nu * c1 / radiance);
   double Tb = (Tbb - beta[var_index]) / alpha[var_index];
   return Tb - 273.15;
}

/** Inverse calculation. Spectral radiance from temperature
* in degree centigrade
 * 
* @param index of the variable to be converted
* @param brightness temperature in [degree C]
* @return radiance value for the given channel
*/
double spectral_radiance(const size_t var_index, const double &temperature) {
   double nu = wavenum[var_index];
   double Tbb = (temperature + 273.15) * alpha[var_index] + beta[var_index];
   return nu * nu * nu * c1 / (exp(c2 * nu / Tbb) - 1);
}

#pragma mark -
#pragma mark MAIN

int main(int argc, char** argv)
{
    // Declare the supported options.
    namespace po = boost::program_options;
    
    po::positional_options_description p;
    p.add("value", -1);
    
    po::options_description desc("Converts spectral radiance to equivalent brightness temperature or vice versa.");
    desc.add_options()
            ("temperature,t", "radiance to temperature")
            ("radiance,r", "temperature to radiance")
            ("value", po::value<double>(), "value to convert")
            ("variable,v", po::value<string>(), "One of msevi_l15_ir_039, msevi_l15_ir_087, msevi_l15_ir_097, msevi_l15_ir_108, msevi_l15_ir_120, msevi_l15_ir_134, msevi_l15_wv_062,msevi_l15_wv_073")
            ;
    
    if (argc < 2) {
        cout << desc << "\n";
        exit(EXIT_SUCCESS);
    }
    
    po::variables_map vm;
    try
    {
        po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
        po::notify(vm);
    } catch (std::exception &e) {
        cerr << "FATAL:parsing command line caused exception: " << e.what()
                << ": check meanie3D-satconv --help for command line options"
                << endl;
        exit(EXIT_FAILURE);
    }

    // Evaluate user input
    double value;
    int variableIndex;
    ConversionType type;
    try {
        parse_commmandline(vm,variableIndex,value,type);
    } catch (const std::exception &e) {
        cerr << "ERROR:exception " << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    switch (type) {
        case RadicanceToTemperature:
            cout << brightness_temperature(variableIndex,value) << endl;
            break;
        case TemperatureToRadiance:
            cout << spectral_radiance(variableIndex,value) << endl;
            break;
    }
    
    return EXIT_SUCCESS;
};
