#include "imt_analysis.h"
#include "onestagepdf_prime.h"
#include "emgpdf.h"
#include "gp_max.h"
#include <math.h>

#ifdef _ENABLE_ONESTAGEPRIME
void
onestagepdf_prime_fixed(const double t[], int t_size, double m, double s,
			double Y[])
{
    double a;
    int k;
    double y;
    double b_y;

    /*  find the value of the first derivitive of the inverse gaussian */
    /*  distribution with parameters m and s at each point in t */
    /*  and return a list of these values corresponding to each point t */
    /* 'onestagepdf_prime:6' Y=exp((-(1-m*t).^2)./(2*s^2*t)); */
    a = 2.0 * (s * s);
    for (k = 0; k < t_size; k++) {
	y = 1.0 - m * t[k];
	Y[k] = exp(-(y * y) / (a * t[k]));
    }

    /* 'onestagepdf_prime:7' Y=(-3/2*t.^(-5/2)+(t.^(-2)/(2*s^2)-m^2/(2*s^2)).*t.^(-3/2))*(1/(s*(2*pi)^.5)).*Y; */
    a = 2.0 * (s * s);
    y = m * m / (2.0 * (s * s));
    b_y = 1.0 / (s * 2.5066282746310002);
    for (k = 0; k < t_size; k++) {
	Y[k] *=
	    (-1.5 * rt_powd_snf(t[k], -2.5) +
	     (rt_powd_snf(t[k], -2.0) / a - y) * rt_powd_snf(t[k],
							     -1.5)) * b_y;
    }
}
#endif

/* End of code generation (onestagepdf_prime.c) */
