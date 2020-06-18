/* 
 * File:   numericalrecipes.h
 * Author: simon
 *
 * Created on October 2, 2014, 1:16 PM
 */

#ifndef NUMERICALRECIPES_H
#define    NUMERICALRECIPES_H

#include <math.h>

#if defined(__cplusplus)
extern "C" {
#endif

float betai(float a, float b, float x);

float betacf(float a, float b, float x);

void crank(unsigned long n, float w[], float *s);

float erfcc(float x);

float gammln(float xx);

void kendl1(float data1[], float data2[],
            unsigned long n, float *tau,
            float *z, float *prob);

void sort2(unsigned long n, float arr[], float brr[]);

void spear(float data1[], float data2[], unsigned long n,
           float *d, float *zd, float *probd, float *rs, float *probrs);

#if defined(__cplusplus)
}
#endif


#endif	

