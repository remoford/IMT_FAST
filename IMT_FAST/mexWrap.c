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

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"



void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{


	/*
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

      mexPrintf ("\n");
    }
	*/

  if (nrhs != 9) {
	  mexPrintf("OMG THE NUMBER OF ARGUMENTS WASNT 8!!!!!\n");
	  return;
  }


  mxClassID zeroth_arg_type = mxGetClassID(prhs[0]);

  if (zeroth_arg_type != mxCHAR_CLASS) {
	  mexPrintf("OMG THE FIRST ARGUMENT WAS NOT A STRING!!!!\n");
	  return;

  }


  // function [P, h, flag, E, C] = convolv_2invG_Dirac_optipon(t, m1, s1, m2, s2, bin, h0, EType)

  if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG THE t ARGUMENT IS NOT OF TYPE DOUBLE!!!!\n");
	  return;
  }

  

  if (mxGetNumberOfElements(prhs[1]) <= 0) {
	  mexPrintf("OMG t needs to have 1 or more elements!!!!\n");
	  return;
  }

  

  if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG THE m1 ARGUMENT IS NOT OF TYPE DOUBLE!!!!\n");
	  return;
  }

  if (mxGetNumberOfElements(prhs[2]) != 1) {
	  mexPrintf("OMG m1 needs to have exactly 1 elements!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG THE s1 ARGUMENT IS NOT OF TYPE DOUBLE!!!!\n");
	  return;
  }

  if (mxGetNumberOfElements(prhs[3]) != 1) {
	  mexPrintf("OMG s1 needs to have exactly 1 elements!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG THE m2 ARGUMENT IS NOT OF TYPE DOUBLE!!!!\n");
	  return;
  }

  if (mxGetNumberOfElements(prhs[4]) != 1) {
	  mexPrintf("OMG m2 needs to have exactly 1 elements!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[5]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG THE s2 ARGUMENT IS NOT OF TYPE DOUBLE!!!!\n");
	  return;
  }

  if (mxGetNumberOfElements(prhs[5]) != 1) {
	  mexPrintf("OMG s2 needs to have exactly 1 elements!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[6]) != mxCHAR_CLASS) {
	  mexPrintf("OMG bin needs to be a char array!!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[7]) != mxDOUBLE_CLASS) {
	  mexPrintf("OMG h0 needs to be of type double!!!!\n");
	  return;
  }

  if (mxGetNumberOfElements(prhs[7]) != 1) {
	  mexPrintf("OMG h0 needs to heave exactly 1 elements!!!!\n");
	  return;
  }

  if (mxGetClassID(prhs[8]) != mxCHAR_CLASS) {
	  mexPrintf("OMG EType needs to be a char array\n");
	  return;
  }


 

  if (nlhs != 1) {
	  mexPrintf("OMG THE NUMBER OF OUTPUTS WASNT 1!!!!!\n");
	  return;
  }

  //if (mxGetClassID(plhs[0]) != mxDOUBLE_CLASS) {
//	  mexPrintf("OMG output needs to be a double type\n");
	//  return;
  //}

  double m1 = mxGetScalar(prhs[2]);
  double s1 = mxGetScalar(prhs[3]);
  double m2 = mxGetScalar(prhs[4]);
  double s2 = mxGetScalar(prhs[5]);

  /* Starting point */
  gsl_vector * v;
  v = gsl_vector_alloc(4);
  gsl_vector_set(v, 0, m1);
  gsl_vector_set(v, 1, s1);
  gsl_vector_set(v, 2, m2);
  gsl_vector_set(v, 3, s2);

  /*
  configStruct config = *(configStruct *)params;

	distType * data = config.data;
	int data_size = config.data_size;
  */

  configStruct config;

  distType * data;
  data = (distType *) mxGetPr(prhs[1]);

  int data_size = mxGetNumberOfElements(prhs[1]);

  config.data = data;
  config.data_size = data_size;



  double ll = convolv_2invG_adapt_nov_loglikelihood(v, (void *) &config);

  plhs[0] = mxCreateDoubleScalar(ll);


  return;


}



#endif