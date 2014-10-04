#ifndef M3D_CONSOLESPINNER_H
#define M3D_CONSOLESPINNER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <iostream>
#include <boost/thread.hpp>

namespace m3D { namespace utils { 

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
        ConsoleSpinner( const ConsoleSpinner &other );

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

        unsigned long                   m_delay_millis;

        vector<char>               m_spinner_chars;

        vector<char>::const_iterator    m_current_position;

        boost::thread                   *m_thread;

        bool                            m_first_char;


#pragma mark - 
#pragma mark Operators

    public:

        /** Used to advance the spinner 
         */
        void
        operator ++();

        /** Used for detaching a thread.
         */
        void
        operator () ();
    };
}}

#endif
