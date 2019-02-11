#include "imt_analysis.h"
#include "twostage.h"
#include "emgpdf.h"
#include "onestagepdf_lag.h"
#include "onestage.h"
#include "tailmass.h"
#include "gsl/gsl_multimin.h"
#include "math.h"
#include "binned_conv.h"
#include "nn_conv.h"
#include "main.h"
#include "loglikelihood.h"
#include "stdio.h"
#include "gsl/gsl_statistics_double.h"
#include "time.h"
#include "utility.h"
#include "gsl/gsl_statistics.h"
#include "mex.h"
#include "matrix.h"


void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  mexPrintf ("Hello World!\n");


  mexPrintf ("nrhs = %d\n", nrhs);
  for (int i = 0; i < nrhs; i++)
    {
      mxClassID classID = mxGetClassID (prhs[i]);
      mexPrintf ("prhs[%d].type = ", i);
      switch (classID)
	{
	case mxINT8_CLASS:
	  mexPrintf ("INT8");
	  break;
	case mxUINT8_CLASS:
	  mexPrintf ("UINT8");
	  break;
	case mxINT16_CLASS:
	  mexPrintf ("INT16");
	  break;
	case mxUINT16_CLASS:
	  mexPrintf ("UINT16");
	  break;
	case mxINT32_CLASS:
	  mexPrintf ("INT32");
	  break;
	case mxUINT32_CLASS:
	  mexPrintf ("UINT32");
	  break;
	case mxINT64_CLASS:
	  mexPrintf ("INT64");
	  break;
	case mxUINT64_CLASS:
	  mexPrintf ("UINT64");
	  break;
	case mxSINGLE_CLASS:
	  mexPrintf ("SINGLE_CLASS");
	  break;
	case mxDOUBLE_CLASS:
	  mexPrintf ("DOUBLE_CLASS");
	  break;
	case mxCHAR_CLASS:
	  mexPrintf ("CHAR_CLASS");
	  break;
	case mxCELL_CLASS:
	  mexPrintf ("CELL_CLASS");
	  break;
	case mxSTRUCT_CLASS:
	  mexPrintf ("STRUCT_CLASS");
	  break;
	case mxUNKNOWN_CLASS:
	  mexPrintf ("UNKNOWN_CLASS");
	  break;
	case mxLOGICAL_CLASS:
	  mexPrintf ("LOGICAL_CLASS");
	  break;
	case mxVOID_CLASS:
	  mexPrintf ("VOID_CLASS");
	  break;
	case mxFUNCTION_CLASS:
	  mexPrintf ("FUNCTION_CLASS");
	  break;
	default:
	  mexPrintf ("UNRECOGNIZED CLASS!!!!");
	}

      mwSize arraySize = mxGetNumberOfElements(prhs[i]);
      mexPrintf(" arraySize = %d", arraySize);

      const char * class_name_str = mxGetClassName(prhs[i]);

      mexPrintf("class name str = %s", class_name_str);

      mexPrintf ("\n");
    }

  return;

  /*
     int ARGSIZE = 1024;
     int argc = 3;
     char **argv = (char **) MALLOC (sizeof (char *));
     for (int i = 0; i < argc; i++)
     {
     argv[i] = (char *) MALLOC (sizeof (char) * ARGSIZE);

     }

     argv[0] = "mexWrap";
     argv[1] = "analysis";
     argv[2] = "onestage";

     main (argc, argv);
   */
}
