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
