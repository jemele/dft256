#include<math.h>
#include "dft.h"
#include"coefficients256.h"

#define USE_LUT

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE])
{
  //Write your code here

	DTYPE pi = 3.141592653589;
    DTYPE outreal[SIZE];
    DTYPE outimag[SIZE];
	int n = SIZE;
    int k;
    for (k = 0; k < n; k++) {  /* For each output element */
        DTYPE sumreal = 0;
        DTYPE sumimag = 0;
        int t;
        for (t = 0; t < n; t++) {  /* For each input element */
#ifdef USE_LUT
            const int angle = (t*k)%n;
            const DTYPE cos_angle = cos_coefficients_table[angle];
            const DTYPE sin_angle = -sin_coefficients_table[angle];
#else
            const DTYPE angle = 0*2 * M_PI * DTYPE(t) * DTYPE(k) / DTYPE(n);
            const DTYPE cos_angle = cos(angle);
            const DTYPE sin_angle = sin(angle);
#endif
            sumreal +=  real_sample[t] * cos_angle + imag_sample[t] * sin_angle;
            sumimag += -real_sample[t] * sin_angle + imag_sample[t] * cos_angle;
        }


        outreal[k] = sumreal;
        outimag[k] = sumimag;
    }

	for(int k = 0; k < SIZE; k++) {
		real_sample[k] = outreal[k];
		imag_sample[k] = outimag[k];
    }
}
