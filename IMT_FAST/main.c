#include "imt_analysis.h"
#include "main.h"
#include "stdio.h"

static void main_IMT_analysis_April2017(const char *);

static void main_IMT_analysis_April2017(const char *model, char * data_filename)
{

	//char filename[1000] = "C:\\Users\\remoford\\source\\repos\\IMT_FAST\\IMT_FAST\\data\\erlotinib_data.txt";
	//char data_filename[1000] = "";
    IMT_analysis_April2017(model, data_filename);
}

int main(int argc, const char *const argv[])
{
    (void) argc;
    (void) argv;

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