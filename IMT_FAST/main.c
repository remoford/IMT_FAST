#include "imt_analysis.h"
#include "main.h"
#include "stdio.h"
#include <math.h>       /* math_errhandling */
#include <errno.h>      /* errno, EDOM */
#include <fenv.h>       /* feclearexcept, fetestexcept, FE_ALL_EXCEPT, FE_INVALID */
#include "float.h"


static void main_IMT_analysis_April2017(const char *, const char * data_filename);

static void main_IMT_analysis_April2017(const char *model, const char * data_filename)
{

	//char filename[1000] = "C:\\Users\\remoford\\source\\repos\\IMT_FAST\\IMT_FAST\\data\\erlotinib_data.txt";
	//char data_filename[1000] = "";
    IMT_analysis_April2017(model, (char *) data_filename);
}

int main(int argc, const char *const argv[])
{
    (void) argc;
    (void) argv;

	/*
	errno = 0;
	if (math_errhandling & MATH_ERREXCEPT) feclearexcept(FE_ALL_EXCEPT);

	printf("Error handling: %d\n", math_errhandling);

	sqrt(-1);
	if (math_errhandling & MATH_ERRNO) {
		if (errno == EDOM) printf("errno set to EDOM\n");
	}
	if (math_errhandling  &MATH_ERREXCEPT) {
		if (fetestexcept(FE_INVALID)) printf("FE_INVALID raised\n");
	}
	*/

	printf("FLT_MIN=%g\n", FLT_MIN);
	printf("FLT_MAX=%g\n", FLT_MAX);
	printf("FLT_MIN_EXP=%g\n", FLT_MIN_EXP);
	printf("FLT_MAX_EXP=%g\n", FLT_MAX_EXP);

	printf("DBL_MIN=%g\n", DBL_MIN);
	printf("DBL_MAX=%g\n", DBL_MAX);
	printf("DBL_MIN_EXP=%g\n", DBL_MIN_EXP);
	printf("DBL_MAX_EXP=%g\n", DBL_MAX_EXP);

	printf("LDBL_MIN=%g\n", LDBL_MIN);
	printf("LDBL_MAX=%g\n", LDBL_MAX);
	printf("LDBL_MIN_EXP=%g\n", LDBL_MIN_EXP);
	printf("LDBL_MAX_EXP=%g\n", LDBL_MAX_EXP);

	/*
	printf("max limit for pow(x,3) : %g\n", pow(DBL_MAX, (1.0 / 3.0)));
	
	printf("second pow limit : %g\n", pow(DBL_MAX, 2.0 / 3.0) / pow(2 * 3.14159, 1.0 / 3.0));
	*/

    /* Invoke the entry-point functions.
       You can call entry-point functions multiple times. */
    if (argc == 2) {
		main_IMT_analysis_April2017(argv[1], "");
	} else if(argc == 3) {
		main_IMT_analysis_April2017(argv[1], argv[2]);
	} else {
		printf("invalid number of arguments, please specify model\n");
    }

    return 0;
}