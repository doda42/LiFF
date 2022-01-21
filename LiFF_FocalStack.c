// LiFF_FocalStack.c 
// 
// mex entry point for C implementation of light field focal stack
// See LiFF_FocalStack.m
// 
// Part of LiFF Light Field Feature Toolbox
// Copyright (c) 2019 Donald G. Dansereau

#include <mex.h>

#include <vl/mathop.h>
#include <vl/sift.h>

#include <math.h>
#include <assert.h>

#include "FocalStack.h"

//---
void mexFunction(int nout, mxArray *out[], int nin, const mxArray *in[])
{
	enum{ IN_LF=0, IN_SLOPEVEC };
	enum{ OUT_FOCSTACK=0 };

	vl_sift_pix const* pLF;
	const mwSize* pLFSize;
	
	double const* pSlopeVec;
	const mwSize* pSlopeVecSize;
	int NumSlopes;

	float* pOutStack = NULL;
	mwSize OutDims[3];

	//--Check args--
	if (nin < 2) {
		mexErrMsgTxt("Two argument required.");
	} else if (nout > 2) {
		mexErrMsgTxt("Too many output arguments.");
	}

	if (mxGetNumberOfDimensions (in[IN_LF]) != 4              ||
		mxGetClassID            (in[IN_LF]) != mxSINGLE_CLASS  ) {
		mexErrMsgTxt("LF must be a monochrome light field of class SINGLE");
	}

	if (mxGetNumberOfDimensions (in[IN_SLOPEVEC]) != 2             ||
		mxGetClassID            (in[IN_SLOPEVEC]) != mxDOUBLE_CLASS) {
		mexErrMsgTxt("SlopeVec must be a vector of DOUBLE");
	}

	//--Grab args--
	pLF = (vl_sift_pix*) mxGetData( in[IN_LF] );
	pLFSize = mxGetDimensions( in[IN_LF] );

	pSlopeVec = (double*) mxGetData( in[IN_SLOPEVEC] );		
	pSlopeVecSize = mxGetDimensions( in[IN_SLOPEVEC] );
	NumSlopes = pSlopeVecSize[1];

	//--Create output array--
	OutDims[0] = NumSlopes;
	OutDims[1] = pLFSize[2];
	OutDims[2] = pLFSize[3];
	out[OUT_FOCSTACK] = mxCreateNumericArray( 3, OutDims, mxSINGLE_CLASS, mxREAL );
	if( out[OUT_FOCSTACK] == NULL )
    	mexErrMsgTxt("Could not create output mxArray.\n");   
    pOutStack = (float*)mxGetData(out[OUT_FOCSTACK]); 

	//--Compute focal stack--
	DoFocalStack( pLF,pLFSize, pSlopeVec,NumSlopes, pOutStack,OutDims );
}
