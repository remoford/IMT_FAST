#include "float.h"
#include "main.h"

void readfile (const distType z[], const distType y[], distType C[], double h, int size_xyz)
{
    for (int i = 0; i < 2*size_xyz; i++) {
		C[i] = 0;
    }

    for (int i = 0; i < size_xyz; i++) {
		for (int j = 0; j < size_xyz; j++) {
			C[i + j] += z[i] * y[j] * h;
		}
    }
}
