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

#ifndef M3D_FILEUTIL_H
#define M3D_FILEUTIL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/filesystem.hpp>

namespace m3D {
    namespace utils {

        /** Finds the part of the filename in the given directory, that
         * is the same for all files. 
         * @param directory
         * @param extension to filter files by (defaults to '*')
         */
        std::string common_component(const std::string &directory,
                                     const std::string extension = "*") {
            namespace fs = boost::filesystem;

            std::string ext = extension;

            if (!boost::algorithm::starts_with(extension, ".")) {
                ext = "." + extension;
            }

            fs::path source_path(directory);

            std::string common;

            if (fs::is_directory(source_path)) {
                fs::directory_iterator dir_iter(source_path);
                fs::directory_iterator end;

                std::vector<std::string> files;

                while (dir_iter != end) {
                    fs::path f = dir_iter->path();

                    if (extension == "*" || boost::algorithm::ends_with(f.filename().generic_string(), ext)) {
                        files.push_back(f.filename().generic_string());
                    }

                    dir_iter++;
                }

                if (!files.empty()) {

                    if (files.size() == 1) {
                        common = files.at(0);
                    } else {
                        size_t token_len = 1;
                        bool have_common = true;
                        while (have_common) {
                            std::string token = files.at(0).substr(0, token_len);
                            size_t num_matching = 0;
                            for (size_t i = 0; i < files.size(); i++) {
                                std::string file = files.at(i);
                                if (boost::algorithm::starts_with(file, token)) {
                                    num_matching++;
                                }
                            }
                            if (num_matching == files.size()) {
                                common = token;
                                token_len++;
                            } else {
                                have_common = false;
                            }
                        }
                    }
                }
            }

            return common;
        }
    }
}

#endif
