
/* Include files */
#include "IMT_analysis_April2017.h"
#include "gp_max.h"
#include "onestagepdf2.h"
#include "emgpdf.h"
#include <gsl/gsl_poly.h>

/* Function Definitions */

double gp_max_fixed(double m, double s)
{
    double coeff[5];
    coeff[4] = pow(m, 4.0) / (4.0 * pow(s, 4.0));
    coeff[3] = 3.0 * (m * m) / (2.0 * (s * s));
    coeff[2] = 3.75 - m * m / (2.0 * pow(s, 4.0));
    coeff[1] = -5.0 / (2.0 * (s * s));
    coeff[0] = 1.0 / (4.0 * pow(s, 4.0));

    double foundRoots[10];

    gsl_poly_complex_workspace *workspace =
	gsl_poly_complex_workspace_alloc(5);
    gsl_poly_complex_solve(coeff, 5, workspace, foundRoots);
    gsl_poly_complex_workspace_free(workspace);

    int numRealRoots = 0;
    double realRoots[10];
    double re;
    double im;
    for (int i = 0; i < 4; i++) {
	re = foundRoots[2 * i];
	im = foundRoots[2 * i + 1];

	if (im <= 0.000000000000001) {
	    if (re >= 0) {
		realRoots[numRealRoots] = re;
		numRealRoots++;
	    }
	}
    }
    double y[10];

    waldpdf(realRoots, m, s, y, numRealRoots);

    double largest = 0;
    for (int i = 0; i < numRealRoots; i++) {
	if (y[i] > largest) {
	    largest = y[i];
	}
    }

    return largest;
}
