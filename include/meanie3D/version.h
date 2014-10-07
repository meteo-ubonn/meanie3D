/* 
 * File:   version.h
 * Author: simon
 *
 * Created on October 7, 2014, 4:01 PM
 */

#ifndef VERSION_H
#define	VERSION_H

#include <string>
#include <iostream>

namespace m3D {
    
    static const std::string VERSION = "1.5.0";
    
    void print_version() {
        std::cout << VERSION << std::endl;
    }
}

#endif	/* VERSION_H */

