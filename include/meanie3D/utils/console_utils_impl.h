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

#ifndef M3D_CONSOLESPINNER_IMPL_H
#define M3D_CONSOLESPINNER_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <sys/time.h>
#include <iostream>
#include <unistd.h>
#include <stdio.h>

#include "console_utils.h"

namespace m3D { namespace utils { 

    ConsoleSpinner::ConsoleSpinner(unsigned long delay_millis)
    : m_delay_millis(delay_millis), m_thread(NULL), m_first_char(true)
    {
        // Simple default set

        m_spinner_chars.push_back('|');
        m_spinner_chars.push_back('/');
        m_spinner_chars.push_back('-');
        m_spinner_chars.push_back('\\');
        m_spinner_chars.push_back('|');
        m_spinner_chars.push_back('/');
        m_spinner_chars.push_back('-');
        m_spinner_chars.push_back('\\');

        m_current_position = m_spinner_chars.begin();
    }


    ConsoleSpinner::ConsoleSpinner(const vector<char> &spinner_chars, unsigned long delay_millis)
    : m_delay_millis(delay_millis),m_spinner_chars(spinner_chars), m_thread(NULL), m_first_char(true)
    {
        m_current_position = m_spinner_chars.begin();
    }

    ConsoleSpinner::ConsoleSpinner( const ConsoleSpinner &o )
    : m_delay_millis( o.m_delay_millis ), m_spinner_chars( o.m_spinner_chars ), m_thread(NULL), m_first_char(true)
    {
        m_current_position = m_spinner_chars.begin();
    }

    ConsoleSpinner::~ConsoleSpinner()
    {
        stop();

        if ( m_thread )
        {
            delete m_thread;
        }
    }

    void
    ConsoleSpinner::start()
    {
        m_thread = new boost::thread( ConsoleSpinner( m_spinner_chars, m_delay_millis ) );

        m_thread->detach();
    }

    void
    ConsoleSpinner::stop()
    {
        if ( m_thread )
        {
            m_thread->interrupt();
        }
    }

    void
    ConsoleSpinner::operator ++()
    {
        if ( m_first_char )
        {
            cout << " ";

            m_first_char = false;
        }

        cout << "\b" << (*m_current_position++);

        m_current_position++;

        if ( m_current_position == m_spinner_chars.end() )
        {
            m_current_position = m_spinner_chars.begin();
        }
    }

    void
    ConsoleSpinner::operator () ()
    {
        while (true)
        {
            this->operator++();

            usleep( m_delay_millis );
        }
    }
}}

#endif
