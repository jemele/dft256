#include<math.h>
#include "dft.h"

static const DTYPE cos_coefficients_table[]={1.000000,0.999699,0.998795,0.997290,0.995185,0.992480,0.989177,0.985278,0.980785,0.975702,0.970031,0.963776,0.956940,0.949528,0.941544,0.932993,0.923880,0.914210,0.903989,0.893224,0.881921,0.870087,0.857729,0.844854,0.831470,0.817585,0.803208,0.788346,0.773010,0.757209,0.740951,0.724247,0.707107,0.689541,0.671559,0.653173,0.634393,0.615232,0.595699,0.575808,0.555570,0.534998,0.514103,0.492898,0.471397,0.449611,0.427555,0.405241,0.382683,0.359895,0.336890,0.313682,0.290285,0.266713,0.242980,0.219101,0.195090,0.170962,0.146730,0.122411,0.098017,0.073565,0.049068,0.024541,0.000000,-0.024541,-0.049068,-0.073565,-0.098017,-0.122411,-0.146730,-0.170962,-0.195090,-0.219101,-0.242980,-0.266713,-0.290285,-0.313682,-0.336890,-0.359895,-0.382683,-0.405241,-0.427555,-0.449611,-0.471397,-0.492898,-0.514103,-0.534998,-0.555570,-0.575808,-0.595699,-0.615232,-0.634393,-0.653173,-0.671559,-0.689541,-0.707107,-0.724247,-0.740951,-0.757209,-0.773010,-0.788346,-0.803208,-0.817585,-0.831470,-0.844854,-0.857729,-0.870087,-0.881921,-0.893224,-0.903989,-0.914210,-0.923880,-0.932993,-0.941544,-0.949528,-0.956940,-0.963776,-0.970031,-0.975702,-0.980785,-0.985278,-0.989177,-0.992480,-0.995185,-0.997290,-0.998795,-0.999699,-1.000000,-0.999699,-0.998795,-0.997290,-0.995185,-0.992480,-0.989177,-0.985278,-0.980785,-0.975702,-0.970031,-0.963776,-0.956940,-0.949528,-0.941544,-0.932993,-0.923880,-0.914210,-0.903989,-0.893224,-0.881921,-0.870087,-0.857729,-0.844854,-0.831470,-0.817585,-0.803208,-0.788346,-0.773010,-0.757209,-0.740951,-0.724247,-0.707107,-0.689541,-0.671559,-0.653173,-0.634393,-0.615232,-0.595699,-0.575808,-0.555570,-0.534998,-0.514103,-0.492898,-0.471397,-0.449611,-0.427555,-0.405241,-0.382683,-0.359895,-0.336890,-0.313682,-0.290285,-0.266713,-0.242980,-0.219101,-0.195090,-0.170962,-0.146730,-0.122411,-0.098017,-0.073565,-0.049068,-0.024541,-0.000000,0.024541,0.049068,0.073565,0.098017,0.122411,0.146730,0.170962,0.195090,0.219101,0.242980,0.266713,0.290285,0.313682,0.336890,0.359895,0.382683,0.405241,0.427555,0.449611,0.471397,0.492898,0.514103,0.534998,0.555570,0.575808,0.595699,0.615232,0.634393,0.653173,0.671559,0.689541,0.707107,0.724247,0.740951,0.757209,0.773010,0.788346,0.803208,0.817585,0.831470,0.844854,0.857729,0.870087,0.881921,0.893224,0.903989,0.914210,0.923880,0.932993,0.941544,0.949528,0.956940,0.963776,0.970031,0.975702,0.980785,0.985278,0.989177,0.992480,0.995185,0.997290,0.998795,0.999699};
static const DTYPE sin_coefficients_table[]={0.000000,-0.024541,-0.049068,-0.073565,-0.098017,-0.122411,-0.146730,-0.170962,-0.195090,-0.219101,-0.242980,-0.266713,-0.290285,-0.313682,-0.336890,-0.359895,-0.382683,-0.405241,-0.427555,-0.449611,-0.471397,-0.492898,-0.514103,-0.534998,-0.555570,-0.575808,-0.595699,-0.615232,-0.634393,-0.653173,-0.671559,-0.689541,-0.707107,-0.724247,-0.740951,-0.757209,-0.773010,-0.788346,-0.803208,-0.817585,-0.831470,-0.844854,-0.857729,-0.870087,-0.881921,-0.893224,-0.903989,-0.914210,-0.923880,-0.932993,-0.941544,-0.949528,-0.956940,-0.963776,-0.970031,-0.975702,-0.980785,-0.985278,-0.989177,-0.992480,-0.995185,-0.997290,-0.998795,-0.999699,-1.000000,-0.999699,-0.998795,-0.997290,-0.995185,-0.992480,-0.989177,-0.985278,-0.980785,-0.975702,-0.970031,-0.963776,-0.956940,-0.949528,-0.941544,-0.932993,-0.923880,-0.914210,-0.903989,-0.893224,-0.881921,-0.870087,-0.857729,-0.844854,-0.831470,-0.817585,-0.803208,-0.788346,-0.773010,-0.757209,-0.740951,-0.724247,-0.707107,-0.689541,-0.671559,-0.653173,-0.634393,-0.615232,-0.595699,-0.575808,-0.555570,-0.534998,-0.514103,-0.492898,-0.471397,-0.449611,-0.427555,-0.405241,-0.382683,-0.359895,-0.336890,-0.313682,-0.290285,-0.266713,-0.242980,-0.219101,-0.195090,-0.170962,-0.146730,-0.122411,-0.098017,-0.073565,-0.049068,-0.024541,-0.000000,0.024541,0.049068,0.073565,0.098017,0.122411,0.146730,0.170962,0.195090,0.219101,0.242980,0.266713,0.290285,0.313682,0.336890,0.359895,0.382683,0.405241,0.427555,0.449611,0.471397,0.492898,0.514103,0.534998,0.555570,0.575808,0.595699,0.615232,0.634393,0.653173,0.671559,0.689541,0.707107,0.724247,0.740951,0.757209,0.773010,0.788346,0.803208,0.817585,0.831470,0.844854,0.857729,0.870087,0.881921,0.893224,0.903989,0.914210,0.923880,0.932993,0.941544,0.949528,0.956940,0.963776,0.970031,0.975702,0.980785,0.985278,0.989177,0.992480,0.995185,0.997290,0.998795,0.999699,1.000000,0.999699,0.998795,0.997290,0.995185,0.992480,0.989177,0.985278,0.980785,0.975702,0.970031,0.963776,0.956940,0.949528,0.941544,0.932993,0.923880,0.914210,0.903989,0.893224,0.881921,0.870087,0.857729,0.844854,0.831470,0.817585,0.803208,0.788346,0.773010,0.757209,0.740951,0.724247,0.707107,0.689541,0.671559,0.653173,0.634393,0.615232,0.595699,0.575808,0.555570,0.534998,0.514103,0.492898,0.471397,0.449611,0.427555,0.405241,0.382683,0.359895,0.336890,0.313682,0.290285,0.266713,0.242980,0.219101,0.195090,0.170962,0.146730,0.122411,0.098017,0.073565,0.049068,0.024541};

#define USE_LUT

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE outreal[SIZE], DTYPE outimag[SIZE])
{
  //Write your code here
	DTYPE pi = 3.141592653589;
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
}
