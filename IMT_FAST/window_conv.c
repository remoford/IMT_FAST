
 /* Include files */
#include "float.h"


void window_conv (const double z[], const double y[], double C[], double h, int size_xyz)
{
    // Assume z is hte concentrated pdf
    //
    int size_conv = 2 * size_xyz;

    for (int i = 0; i < size_conv; i++) {
	C[i] = 0;
    }


    /* Find the highest zero-valued index */
    int firstIdx = 0;
    for (int i = 0; i < size_xyz; i++) {
        if (z[i] > 0)
		break;
        else
            firstIdx = i;

    }




    /* do the lopsided convolution */
    for (int i = firstIdx; i < size_xyz; i++) {
	for (int j = 0; j < size_xyz; j++) {
	    C[i + j] += z[i] * y[j] * h;
	}
    }






}
