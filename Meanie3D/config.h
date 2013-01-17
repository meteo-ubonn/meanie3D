#ifndef CF_ALGORITHMS_CONFIG_H
#define CF_ALGORITHMS_CONFIG_H

/* OS - type */
/* #undef cf */
/* #undef cf */

#if defined(CF_ALGORITHMS_OS_MACOSX) || defined(CF_ALGORITHMS_OS_LINUX)
#define CF_ALGORITHMS_OS_POSIX
#endif

#endif /* header guard */
