#include "imt_analysis.h"
#include "main.h"
#include "stdio.h"

static void main_IMT_analysis_April2017(const char *);

static void main_IMT_analysis_April2017(const char *model)
{
    IMT_analysis_April2017(model);
}

int main(int argc, const char *const argv[])
{
    (void) argc;
    (void) argv;

    /* Invoke the entry-point functions.
       You can call entry-point functions multiple times. */
    if (argc == 2) {
	main_IMT_analysis_April2017(argv[1]);
    } else {
	printf("invalid number of arguments, please specify model\n");
    }

    return 0;
}