/* 
 * File:   parallel.h
 * Author: simon
 *
 * Created on September 15, 2014, 9:52 PM
 */

#ifndef PARALL_H
#define	PARALL_H

#if WITH_TBB
#include <tbb/tbb.h>
#endif

#if WITH_BOOST_THREADS
#include <boost/thread.hpp>
#endif

#if WITH_OPENMP
#include <omp.h>
#endif

#endif	/* PARALLEL_H */

