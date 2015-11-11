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

#ifndef M3D_CONSOLESPINNER_H
#define M3D_CONSOLESPINNER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <iostream>
#include <boost/thread.hpp>

namespace m3D {
    namespace utils {

        class ConsoleSpinner
        {
        public:

#pragma mark -
#pragma mark Constructors / Destructors

            /** Constructs a spinner.
             * @param delay between spinner updates in milliseconds.
             */
            ConsoleSpinner(unsigned long delay_millis = 250);

            /** Constructs a spinner.
             * @param give the spinner your own set of characters to spin through
             * @param delay between spinner updates in milliseconds.
             */
            ConsoleSpinner(const vector<char> &spinner_chars, unsigned long delay_millis = 250);

            /** Stops the spinner as well
             */
            ~ConsoleSpinner();

            /** Copy constructor
             * @param console spinner
             */
            ConsoleSpinner(const ConsoleSpinner &other);

#pragma mark -
#pragma mark Starting and Stopping

            /** Starts the spinner
             */
            void start();

            /** Stops the spinner 
             */
            void stop();

        private:

#pragma mark -
#pragma mark Private Member Variables

            unsigned long m_delay_millis;

            vector<char> m_spinner_chars;

            vector<char>::const_iterator m_current_position;

            boost::thread *m_thread;

            bool m_first_char;


#pragma mark - 
#pragma mark Operators

        public:

            /** Used to advance the spinner 
             */
            void
            operator++();

            /** Used for detaching a thread.
             */
            void
            operator()();
        };
    }
}

#endif
