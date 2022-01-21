/*
LiFF_ExtractFeatures.c

mex entry point for C implementation of LiFF light field features
See LiFF_ExtractFeatures.m

A full description of the LiFF feature detector and descriptor are available here:
[1] D. G. Dansereau, B. Girod, and G. Wetzstein, “LiFF: Light field features in scale and depth,” 
in Computer Vision and Pattern Recognition (CVPR), 2019. 
Paper and supplemental information are at https://roboticimaging.org/Tools/LiFF/

Part of LiFF Light Field Feature Toolbox
Copyright (c) 2019 Donald G. Dansereau

Based on code from the VLFeat library 
Copyright (C) 2007-12 by Andrea Vedaldi and Brian Fulkerson
Released under the terms of the BSD license, see COPYING.txt.
*/

#include <toolbox/mexutils.h>
#include <vl/mathop.h>
#include <vl/sift.h>

#include <math.h>
#include <assert.h>

#include "FocalStack.h"

//---Fwd decl---
void LiFF_detect(VlSiftFilt *fs[], int num_slopes);
VL_INLINE void L1Root_normalize(vl_sift_pix *begin, vl_sift_pix *end);

/* option codes */
enum
{
	opt_octaves = 0,
	opt_levels,
	opt_first_octave,
	opt_frames,
	opt_edge_thresh,
	opt_peak_thresh,
	opt_norm_thresh,
	opt_magnif,
	opt_window_size,
	opt_orientations,
	opt_float_descriptors,
	opt_first_slope,
	opt_last_slope,
	opt_num_slopes,
	opt_L2_norm,
	opt_verbose
};

/* options */
vlmxOption options[] = {
	{"Octaves", 1, opt_octaves},
	{"Levels", 1, opt_levels},
	{"FirstOctave", 1, opt_first_octave},
	{"Frames", 1, opt_frames},
	{"PeakThresh", 1, opt_peak_thresh},
	{"EdgeThresh", 1, opt_edge_thresh},
	{"NormThresh", 1, opt_norm_thresh},
	{"Magnif", 1, opt_magnif},
	{"WindowSize", 1, opt_window_size},
	{"FirstSlope", 1, opt_first_slope},
	{"LastSlope", 1, opt_last_slope},
	{"NumSlopes", 1, opt_num_slopes},
	{"Orientations", 0, opt_orientations},
	{"FloatDescriptors", 0, opt_float_descriptors},
	{"L2Norm", 0, opt_L2_norm},
	{"Verbose", 0, opt_verbose},
	{0, 0, 0}};

#define CHECK_NEIGHBORS(v, pt, CMP, SGN)      \
	(v CMP## = SGN 0.8 * tp &&                \
			   v CMP * (pt + xo) &&           \
			   v CMP * (pt - xo) &&           \
			   v CMP * (pt + so) &&           \
			   v CMP * (pt - so) &&           \
			   v CMP * (pt + yo) &&           \
			   v CMP * (pt - yo) &&           \
                                              \
			   v CMP * (pt + yo + xo) &&      \
			   v CMP * (pt + yo - xo) &&      \
			   v CMP * (pt - yo + xo) &&      \
			   v CMP * (pt - yo - xo) &&      \
                                              \
			   v CMP * (pt + xo + so) &&      \
			   v CMP * (pt - xo + so) &&      \
			   v CMP * (pt + yo + so) &&      \
			   v CMP * (pt - yo + so) &&      \
			   v CMP * (pt + yo + xo + so) && \
			   v CMP * (pt + yo - xo + so) && \
			   v CMP * (pt - yo + xo + so) && \
			   v CMP * (pt - yo - xo + so) && \
                                              \
			   v CMP * (pt + xo - so) &&      \
			   v CMP * (pt - xo - so) &&      \
			   v CMP * (pt + yo - so) &&      \
			   v CMP * (pt - yo - so) &&      \
			   v CMP * (pt + yo + xo - so) && \
			   v CMP * (pt + yo - xo - so) && \
			   v CMP * (pt - yo + xo - so) && \
			   v CMP * (pt - yo - xo - so))

#define CHECK_NEIGHBORS4(v, pt, ptNext, ptPrev, CMP, SGN)                                \
	(CHECK_NEIGHBORS(v, pt, CMP, SGN) &&                                                 \
	 ((ptNext == pt) || ((v CMP * (ptNext)) && CHECK_NEIGHBORS(v, ptNext, CMP, SGN))) && \
	 ((ptPrev == pt) || ((v CMP * (ptPrev)) && CHECK_NEIGHBORS(v, ptPrev, CMP, SGN))))

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Transpose desriptor
 **
 ** @param dst destination buffer.
 ** @param src source buffer.
 **
 ** The function writes to @a dst the transpose of the SIFT descriptor
 ** @a src. The tranpsose is defined as the descriptor that one
 ** obtains from computing the normal descriptor on the transposed
 ** image.
 **/

VL_INLINE void
transpose_descriptor(vl_sift_pix *dst, vl_sift_pix *src)
{
	int const BO = 8; /* number of orientation bins */
	int const BP = 4; /* number of spatial bins	 */
	int i, j, t;

	for (j = 0; j < BP; ++j)
	{
		int jp = BP - 1 - j;
		for (i = 0; i < BP; ++i)
		{
			int o = BO * i + BP * BO * j;
			int op = BO * i + BP * BO * jp;
			dst[op] = src[o];
			for (t = 1; t < BO; ++t)
				dst[BO - t + op] = src[t + o];
		}
	}
}

// Copy a slice out of the focal stack into an image that can be passed into the existing SIFT code
void LiFF_CopySlice(vl_sift_pix *pDest, int iSlope, const vl_sift_pix *pFocStack, const mwSize *pStackSize)
{
	for (int VIdx = 0; VIdx < pStackSize[1]; ++VIdx)
	{
		for (int UIdx = 0; UIdx < pStackSize[2]; ++UIdx)
		{
			vl_sift_pix CurVal = pFocStack[iSlope + VIdx * pStackSize[0] + UIdx * pStackSize[0] * pStackSize[1]];
			pDest[VIdx + UIdx * pStackSize[1]] = CurVal;
		}
	}
}

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Ordering of tuples by increasing scale
 **
 ** @param a tuple.
 ** @param b tuple.
 **
 ** @return @c a[2] < b[2]
 **/

static int
korder(void const *a, void const *b)
{
	double x = ((double *)a)[2] - ((double *)b)[2];
	if (x < 0)
		return -1;
	if (x > 0)
		return +1;
	return 0;
}

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Check for sorted keypoints
 **
 ** @param keys keypoint list to check
 ** @param nkeys size of the list.
 **
 ** @return 1 if the keypoints are storted.
 **/

vl_bool
check_sorted(double const *keys, vl_size nkeys)
{
	vl_uindex k;
	for (k = 0; k + 1 < nkeys; ++k)
	{
		if (korder(keys, keys + 5) > 0)
		{
			return VL_FALSE;
		}
		keys += 5;
	}
	return VL_TRUE;
}

/** ------------------------------------------------------------------
 ** @brief MEX entry point
 **/

void mexFunction(int nout, mxArray *out[],
				 int nin, const mxArray *in[])
{
	enum
	{
		IN_LF = 0,
		IN_END
	};
	enum
	{
		OUT_FRAMES = 0,
		OUT_DESCRIPTORS,
		OUT_FOCSTACK
	};

	int verbose = 0;
	int L2Norm = 0;
	int opt;
	int next = IN_END;
	mxArray const *optarg;

	int O = -1;
	int S = 3;
	int o_min = 0;

	double first_slope = -1.0;
	double last_slope = 1.0;
	int num_slopes = -1;

	double edge_thresh = -1;
	double peak_thresh = -1;
	double norm_thresh = -1;
	double magnif = -1;
	double window_size = -1;

	mxArray *ikeys_array = 0;
	double *ikeys = 0;
	int nikeys = -1;
	vl_bool force_orientations = 0;
	vl_bool floatDescriptors = 0;

	vl_sift_pix const *pLF = NULL;
	const mwSize *pLFSize = NULL;

	VL_USE_MATLAB_ENV;

	/* -----------------------------------------------------------------
	 *	 Check the arguments
	 * -------------------------------------------------------------- */

	if (nin < 1)
	{
		mexErrMsgTxt("One argument required.");
	}
	else if (nout > 3)
	{
		mexErrMsgTxt("Too many output arguments.");
	}

	if (mxGetNumberOfDimensions(in[IN_LF]) != 4 ||
		mxGetClassID(in[IN_LF]) != mxSINGLE_CLASS)
	{
		mexErrMsgTxt("LF must be a monochrome light field of class SINGLE");
	}

	//--Grab args--
	pLF = (vl_sift_pix *)mxGetData(in[IN_LF]);
	pLFSize = mxGetDimensions(in[IN_LF]);

	while ((opt = vlmxNextOption(in, nin, options, &next, &optarg)) >= 0)
	{
		switch (opt)
		{

		case opt_first_slope:
			if (!vlmxIsPlainScalar(optarg))
			{
				mexErrMsgTxt("'FirstSlope' must be a real.");
			}
			first_slope = *mxGetPr(optarg);
			break;
		case opt_last_slope:
			if (!vlmxIsPlainScalar(optarg))
			{
				mexErrMsgTxt("'LastSlope' must be a real.");
			}
			last_slope = *mxGetPr(optarg);
			break;
		case opt_num_slopes:
			if (!vlmxIsPlainScalar(optarg) || (num_slopes = (int)*mxGetPr(optarg)) < 1)
			{
				mexErrMsgTxt("'NumSlopes' must be a positive integer.");
			}
			break;

		case opt_L2_norm:
			L2Norm = 1;
			break;

		case opt_verbose:
			++verbose;
			break;

		case opt_octaves:
			if (!vlmxIsPlainScalar(optarg) || (O = (int)*mxGetPr(optarg)) < 0)
			{
				mexErrMsgTxt("'Octaves' must be a positive integer.");
			}
			break;

		case opt_levels:
			if (!vlmxIsPlainScalar(optarg) || (S = (int)*mxGetPr(optarg)) < 1)
			{
				mexErrMsgTxt("'Levels' must be a positive integer.");
			}
			break;

		case opt_first_octave:
			if (!vlmxIsPlainScalar(optarg))
			{
				mexErrMsgTxt("'FirstOctave' must be an integer");
			}
			o_min = (int)*mxGetPr(optarg);
			break;

		case opt_edge_thresh:
			if (!vlmxIsPlainScalar(optarg) || (edge_thresh = *mxGetPr(optarg)) < 1)
			{
				mexErrMsgTxt("'EdgeThresh' must be not smaller than 1.");
			}
			break;

		case opt_peak_thresh:
			if (!vlmxIsPlainScalar(optarg) || (peak_thresh = *mxGetPr(optarg)) < 0)
			{
				mexErrMsgTxt("'PeakThresh' must be a non-negative real.");
			}
			break;

		case opt_norm_thresh:
			if (!vlmxIsPlainScalar(optarg) || (norm_thresh = *mxGetPr(optarg)) < 0)
			{
				mexErrMsgTxt("'NormThresh' must be a non-negative real.");
			}
			break;

		case opt_magnif:
			if (!vlmxIsPlainScalar(optarg) || (magnif = *mxGetPr(optarg)) < 0)
			{
				mexErrMsgTxt("'Magnif' must be a non-negative real.");
			}
			break;

		case opt_window_size:
			if (!vlmxIsPlainScalar(optarg) || (window_size = *mxGetPr(optarg)) < 0)
			{
				mexErrMsgTxt("'WindowSize' must be a non-negative real.");
			}
			break;

		case opt_frames:
			if (!vlmxIsMatrix(optarg, 5, -1))
			{
				mexErrMsgTxt("'Frames' must be a 5 x N matrix.");
			}
			ikeys_array = mxDuplicateArray(optarg);
			nikeys = mxGetN(optarg);
			ikeys = mxGetPr(ikeys_array);
			if (!check_sorted(ikeys, nikeys))
			{
				qsort(ikeys, nikeys, 5 * sizeof(double), korder);
			}
			break;

		case opt_orientations:
			force_orientations = 1;
			break;

		case opt_float_descriptors:
			floatDescriptors = 1;
			break;

		default:
			abort();
		}
	}

	//--Defaults-------------
	if (num_slopes <= 0)
	{
		int MaxSpatialSamps = (pLFSize[0] > pLFSize[1]) ? pLFSize[0] : pLFSize[1];
		num_slopes = (MaxSpatialSamps / 2) * 2 + 1; // force odd
	}

	//--Find focal stack-------------
	double *pSlopeVec = NULL;
	pSlopeVec = mxMalloc(num_slopes * sizeof(double));

	double CurSlope = first_slope;
	double SlopeStep = (last_slope - first_slope) / (num_slopes - 1);
	for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
	{
		pSlopeVec[iSlope] = CurSlope;
		CurSlope += SlopeStep;
	}

	vl_sift_pix *pFocStack = NULL;
	mwSize pStackSize[3] = {num_slopes, pLFSize[2], pLFSize[3]};
	mxArray *pFocStackOut = mxCreateNumericArray(3, pStackSize, mxSINGLE_CLASS, mxREAL);
	if (pFocStackOut == NULL)
		mexErrMsgTxt("Could not create output mxArray.\n");
	pFocStack = (float *)mxGetData(pFocStackOut);

	DoFocalStack(pLF, pLFSize, pSlopeVec, num_slopes, pFocStack, pStackSize);

	//--allocate space for a single view-----------------------------
	vl_sift_pix *data = NULL;
	data = mxMalloc(pLFSize[2] * pLFSize[3] * sizeof(vl_sift_pix));

	int M = pLFSize[2];
	int N = pLFSize[3];

	//-----------------------------------------------------------------
	{
		VlSiftFilt *filts[num_slopes];
		vl_bool first;
		double *frames = 0;
		void *descr = 0;
		int nframes = 0, reserved = 0, i, j, q;

		for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
		{
			/* create a filter to process the image */
			filts[iSlope] = vl_sift_new(M, N, O, S, o_min);

			if (peak_thresh >= 0)
				vl_sift_set_peak_thresh(filts[iSlope], peak_thresh);
			if (edge_thresh >= 0)
				vl_sift_set_edge_thresh(filts[iSlope], edge_thresh);
			if (norm_thresh >= 0)
				vl_sift_set_norm_thresh(filts[iSlope], norm_thresh);
			if (magnif >= 0)
				vl_sift_set_magnif(filts[iSlope], magnif);
			if (window_size >= 0)
				vl_sift_set_window_size(filts[iSlope], window_size);
		}

		if (verbose)
		{
			mexPrintf("vl_sift: filter settings:\n");
			mexPrintf("vl_sift:	 octaves			(O)			= %d\n",
					  vl_sift_get_noctaves(filts[0]));
			mexPrintf("vl_sift:	 levels			 (S)			= %d\n",
					  vl_sift_get_nlevels(filts[0]));
			mexPrintf("vl_sift:	 first octave (o_min)	= %d\n",
					  vl_sift_get_octave_first(filts[0]));
			mexPrintf("vl_sift:	 edge thresh					 = %g\n",
					  vl_sift_get_edge_thresh(filts[0]));
			mexPrintf("vl_sift:	 peak thresh					 = %g\n",
					  vl_sift_get_peak_thresh(filts[0]));
			mexPrintf("vl_sift:	 norm thresh					 = %g\n",
					  vl_sift_get_norm_thresh(filts[0]));
			mexPrintf("vl_sift:	 window size					 = %g\n",
					  vl_sift_get_window_size(filts[0]));
			mexPrintf("vl_sift:	 float descriptor			= %d\n",
					  floatDescriptors);

			mexPrintf("vl_sift:	 first slope	= %g\n", first_slope);
			mexPrintf("vl_sift:	 last slope		= %g\n", last_slope);
			mexPrintf("vl_sift:	 num slopes		= %d\n", num_slopes);
			mexPrintf("vl_sift:	 Normalizing descriptors using %s-norm\n", L2Norm ? "L2" : "L1 root");

			mexPrintf((nikeys >= 0) ? "vl_sift: will source frames? yes (%d read)\n" : "vl_sift: will source frames? no\n", nikeys);
			mexPrintf("vl_sift: will force orientations? %s\n",
					  force_orientations ? "yes" : "no");
		}

		/* ...............................................................
		 *	Process each octave
		 * ............................................................ */
		i = 0;
		first = 1;
		while (1)
		{
			int err;
			VlSiftKeypoint const *keys = 0;
			int nkeys = 0;

			if (verbose)
			{
				mexPrintf("vl_sift: processing octave %d\n",
						  vl_sift_get_octave_index(filts[0]));
			}

			/* Calculate the GSS for the next octave .................... */
			for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
			{
				if (first)
				{
					LiFF_CopySlice(data, iSlope, pFocStack, pStackSize);
					err = vl_sift_process_first_octave(filts[iSlope], data);
				}
				else
				{
					err = vl_sift_process_next_octave(filts[iSlope]);
				}
			}
			first = 0;
			if (err)
				break;

			if (verbose > 1)
			{
				mexPrintf("vl_sift: GSS octave %d computed\n",
						  vl_sift_get_octave_index(filts[0]));
			}

			/* Run detector ............................................. */
			if (nikeys < 0)
			{
				LiFF_detect(filts, num_slopes);
				nkeys = -1; // we will load batches of keys, one batch per slope
			}
			else
			{
				nkeys = nikeys; // keys provided as argument
			}

			/* For each keypoint ........................................ */
			int CurSlopeIdx = -1;
			int CurKeyIdx = 0;
			while (1)
			{
				double angles[4];
				int nangles;
				VlSiftKeypoint ik;
				VlSiftKeypoint const *k;

				/* Obtain keypoint orientations ........................... */
				if (nikeys >= 0)
				{
					// Keypoints were provided as an input argument
					if (CurKeyIdx >= nikeys)
						break;

					CurSlopeIdx = (int)ikeys[5 * CurKeyIdx + 4] - 1;
					if (CurSlopeIdx < 0 || CurSlopeIdx >= num_slopes)
						mexErrMsgTxt("Feature's input slope is out of range");
					vl_sift_keypoint_init(filts[CurSlopeIdx], &ik,
										  ikeys[5 * CurKeyIdx + 1] - 1,
										  ikeys[5 * CurKeyIdx + 0] - 1,
										  ikeys[5 * CurKeyIdx + 2]);

					if (ik.o != vl_sift_get_octave_index(filts[CurSlopeIdx]))
					{
						++CurKeyIdx;
						continue; // keypoint belongs to a different octave
					}
					k = &ik;

					/* optionally compute orientations too */
					if (force_orientations)
					{
						nangles = vl_sift_calc_keypoint_orientations(filts[CurSlopeIdx], angles, k);
					}
					else
					{
						angles[0] = VL_PI / 2 - ikeys[5 * CurKeyIdx + 3];
						nangles = 1;
					}
					++CurKeyIdx;
				}
				else
				{
					// grab keys from the next slope
					while (CurKeyIdx >= nkeys && CurSlopeIdx < num_slopes)
					{
						++CurSlopeIdx;
						if (CurSlopeIdx >= num_slopes)
							break;
						keys = vl_sift_get_keypoints(filts[CurSlopeIdx]);
						nkeys = vl_sift_get_nkeypoints(filts[CurSlopeIdx]);
						CurKeyIdx = 0;
					}
					if (CurSlopeIdx >= num_slopes)
						break;

					k = keys + CurKeyIdx;
					++CurKeyIdx;
					nangles = vl_sift_calc_keypoint_orientations(filts[CurSlopeIdx], angles, k);
				}

				/* For each orientation ................................... */
				for (q = 0; q < nangles; ++q)
				{
					vl_sift_pix buf[128];
					vl_sift_pix rbuf[128];

					/* compute descriptor (if necessary) */
					if (nout > 1)
					{
						vl_sift_calc_keypoint_descriptor(filts[CurSlopeIdx], buf, k, angles[q]);
						if (!L2Norm)
						{
#define NBO 8
#define NBP 4
							L1Root_normalize(buf, buf + NBO * NBP * NBP);
						}
						transpose_descriptor(rbuf, buf);
					}

					/* make enough room for all these keypoints and more */
					if (reserved < nframes + 1)
					{
						reserved += 2 * nkeys;
						frames = mxRealloc(frames, 5 * sizeof(double) * reserved);
						if (nout > 1)
						{
							if (!floatDescriptors)
							{
								descr = mxRealloc(descr, 128 * sizeof(vl_uint8) * reserved);
							}
							else
							{
								descr = mxRealloc(descr, 128 * sizeof(float) * reserved);
							}
						}
					}

					/* Save back with MATLAB conventions. Notice tha the input
					 * image was the transpose of the actual image. */
					frames[5 * nframes + 0] = k->y + 1;
					frames[5 * nframes + 1] = k->x + 1;
					frames[5 * nframes + 2] = k->sigma;
					frames[5 * nframes + 3] = VL_PI / 2 - angles[q];
					frames[5 * nframes + 4] = CurSlopeIdx + 1;

					if (nout > 1)
					{
						if (!floatDescriptors)
						{
							for (j = 0; j < 128; ++j)
							{
								float x = 512.0F * rbuf[j];
								x = (x < 255.0F) ? x : 255.0F;
								((vl_uint8 *)descr)[128 * nframes + j] = (vl_uint8)x;
							}
						}
						else
						{
							for (j = 0; j < 128; ++j)
							{
								float x = 512.0F * rbuf[j];
								((float *)descr)[128 * nframes + j] = x;
							}
						}
					}

					++nframes;
				} /* next orientation */
			}	 /* next keypoint */
		}		  /* next octave */

		if (verbose)
		{
			mexPrintf("vl_sift: found %d keypoints\n", nframes);
		}

		/* ...............................................................
		 *	Save back
		 * ............................................................ */

		{
			mwSize dims[2];

			/* create an empty array */
			dims[0] = 0;
			dims[1] = 0;
			out[OUT_FRAMES] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

			/* set array content to be the frames buffer */
			dims[0] = 5;
			dims[1] = nframes;
			mxSetPr(out[OUT_FRAMES], frames);
			mxSetDimensions(out[OUT_FRAMES], dims, 2);

			if (nout > 1)
			{
				/* create an empty array */
				dims[0] = 0;
				dims[1] = 0;
				out[OUT_DESCRIPTORS] = mxCreateNumericArray(2, dims,
															floatDescriptors ? mxSINGLE_CLASS : mxUINT8_CLASS,
															mxREAL);

				/* set array content to be the descriptors buffer */
				dims[0] = 128;
				dims[1] = nframes;
				mxSetData(out[OUT_DESCRIPTORS], descr);
				mxSetDimensions(out[OUT_DESCRIPTORS], dims, 2);
			}

			if (nout > 2)
			{
				out[OUT_FOCSTACK] = pFocStackOut;
			}
		}

		/* cleanup */
		for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
			vl_sift_delete(filts[iSlope]);

		if (ikeys_array)
			mxDestroyArray(ikeys_array);

	} /* end: do job */

	if (data)
		mxFree(data);

	if (pSlopeVec)
		mxFree(pSlopeVec);
}

// Modified sift_detect
void LiFF_detect(VlSiftFilt *fs[], int num_slopes)
{
	VlSiftFilt *f = fs[0];
	int s_min = f->s_min;
	int s_max = f->s_max;
	int w = f->octave_width;
	int h = f->octave_height;
	double te = f->edge_thresh;
	double tp = f->peak_thresh;

	int const xo = 1;	 /* x-stride */
	int const yo = w;	 /* y-stride */
	int const so = w * h; /* s-stride */

	double xper = pow(2.0, f->o_cur);

	int x, y, s, i, ii, jj;
	vl_sift_pix *pt, v;
	VlSiftKeypoint *k;

	/* compute difference of gaussian (DoG) */
	for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
	{
		fs[iSlope]->nkeys = 0; // clear list of keys

		pt = fs[iSlope]->dog;
		for (s = s_min; s <= s_max - 1; ++s)
		{
			vl_sift_pix *src_a = vl_sift_get_octave(fs[iSlope], s);
			vl_sift_pix *src_b = vl_sift_get_octave(fs[iSlope], s + 1);
			vl_sift_pix *end_a = src_a + w * h;
			while (src_a != end_a)
			{
				*pt++ = *src_b++ - *src_a++;
			}
		}
	}

	//--Find local maxima of DoG---------------------------------------------------
	for (int iSlope = 0; iSlope < num_slopes; ++iSlope)
	{
		f = fs[iSlope];
		vl_sift_pix *dog = f->dog;

		/* start from dog [1,1,s_min+1] */
		pt = fs[iSlope]->dog + xo + yo + so;

		vl_sift_pix *ptNext = pt;
		vl_sift_pix *ptPrev = pt;

		if (iSlope > 0)
			ptPrev = fs[iSlope - 1]->dog + xo + yo + so;

		if (iSlope < num_slopes - 1)
			ptNext = fs[iSlope + 1]->dog + xo + yo + so;

		for (s = s_min + 1; s <= s_max - 2; ++s)
		{
			for (y = 1; y < h - 1; ++y)
			{
				for (x = 1; x < w - 1; ++x)
				{
					v = *pt;
					if ((CHECK_NEIGHBORS4(v, pt, ptNext, ptPrev, >, +)) ||
						(CHECK_NEIGHBORS4(v, pt, ptNext, ptPrev, <, -)))
					{

						/* make room for more keypoints */
						if (f->nkeys >= f->keys_res)
						{
							f->keys_res += 500;
							if (f->keys)
							{
								f->keys = vl_realloc(f->keys,
													 f->keys_res *
														 sizeof(VlSiftKeypoint));
							}
							else
							{
								f->keys = vl_malloc(f->keys_res *
													sizeof(VlSiftKeypoint));
							}
						}

						k = f->keys + (f->nkeys++);

						k->ix = x;
						k->iy = y;
						k->is = s;
					}
					pt += 1;
					ptNext += 1;
					ptPrev += 1;
				}
				pt += 2;
				ptNext += 2;
				ptPrev += 2;
			}
			pt += 2 * yo;
			ptNext += 2 * yo;
			ptPrev += 2 * yo;
		}

		/* -----------------------------------------------------------------
		*                                               Refine local maxima
		* -------------------------------------------------------------- */

		/* this pointer is used to write the keypoints back */
		k = f->keys;

		for (i = 0; i < f->nkeys; ++i)
		{

			int x = f->keys[i].ix;
			int y = f->keys[i].iy;
			int s = f->keys[i].is;

			double Dx = 0, Dy = 0, Ds = 0, Dxx = 0, Dyy = 0, Dss = 0, Dxy = 0, Dxs = 0, Dys = 0;
			double A[3 * 3], b[3];

			int dx = 0;
			int dy = 0;

			int iter, i, j;

			for (iter = 0; iter < 5; ++iter)
			{

				x += dx;
				y += dy;

				pt = fs[iSlope]->dog + xo * x + yo * y + so * (s - s_min);

				/** @brief Index GSS @internal */
#define at(dx, dy, ds) (*(pt + (dx)*xo + (dy)*yo + (ds)*so))

				/** @brief Index matrix A @internal */
#define Aat(i, j) (A[(i) + (j)*3])

				/* compute the gradient */
				Dx = 0.5 * (at(+1, 0, 0) - at(-1, 0, 0));
				Dy = 0.5 * (at(0, +1, 0) - at(0, -1, 0));
				Ds = 0.5 * (at(0, 0, +1) - at(0, 0, -1));

				/* compute the Hessian */
				Dxx = (at(+1, 0, 0) + at(-1, 0, 0) - 2.0 * at(0, 0, 0));
				Dyy = (at(0, +1, 0) + at(0, -1, 0) - 2.0 * at(0, 0, 0));
				Dss = (at(0, 0, +1) + at(0, 0, -1) - 2.0 * at(0, 0, 0));

				Dxy = 0.25 * (at(+1, +1, 0) + at(-1, -1, 0) - at(-1, +1, 0) - at(+1, -1, 0));
				Dxs = 0.25 * (at(+1, 0, +1) + at(-1, 0, -1) - at(-1, 0, +1) - at(+1, 0, -1));
				Dys = 0.25 * (at(0, +1, +1) + at(0, -1, -1) - at(0, -1, +1) - at(0, +1, -1));

				/* solve linear system ....................................... */
				Aat(0, 0) = Dxx;
				Aat(1, 1) = Dyy;
				Aat(2, 2) = Dss;
				Aat(0, 1) = Aat(1, 0) = Dxy;
				Aat(0, 2) = Aat(2, 0) = Dxs;
				Aat(1, 2) = Aat(2, 1) = Dys;

				b[0] = -Dx;
				b[1] = -Dy;
				b[2] = -Ds;

				/* Gauss elimination */
				for (j = 0; j < 3; ++j)
				{
					double maxa = 0;
					double maxabsa = 0;
					int maxi = -1;
					double tmp;

					/* look for the maximally stable pivot */
					for (i = j; i < 3; ++i)
					{
						double a = Aat(i, j);
						double absa = vl_abs_d(a);
						if (absa > maxabsa)
						{
							maxa = a;
							maxabsa = absa;
							maxi = i;
						}
					}

					/* if singular give up */
					if (maxabsa < 1e-10f)
					{
						b[0] = 0;
						b[1] = 0;
						b[2] = 0;
						break;
					}

					i = maxi;

					/* swap j-th row with i-th row and normalize j-th row */
					for (jj = j; jj < 3; ++jj)
					{
						tmp = Aat(i, jj);
						Aat(i, jj) = Aat(j, jj);
						Aat(j, jj) = tmp;
						Aat(j, jj) /= maxa;
					}
					tmp = b[j];
					b[j] = b[i];
					b[i] = tmp;
					b[j] /= maxa;

					/* elimination */
					for (ii = j + 1; ii < 3; ++ii)
					{
						double x = Aat(ii, j);
						for (jj = j; jj < 3; ++jj)
						{
							Aat(ii, jj) -= x * Aat(j, jj);
						}
						b[ii] -= x * b[j];
					}
				}

				/* backward substitution */
				for (i = 2; i > 0; --i)
				{
					double x = b[i];
					for (ii = i - 1; ii >= 0; --ii)
					{
						b[ii] -= x * Aat(ii, i);
					}
				}

				/* .......................................................... */
				/* If the translation of the keypoint is big, move the keypoint
				* and re-iterate the computation. Otherwise we are all set.
				*/

				dx = ((b[0] > 0.6 && x < w - 2) ? 1 : 0) + ((b[0] < -0.6 && x > 1) ? -1 : 0);

				dy = ((b[1] > 0.6 && y < h - 2) ? 1 : 0) + ((b[1] < -0.6 && y > 1) ? -1 : 0);

				if (dx == 0 && dy == 0)
					break;
			}

			/* check threshold and other conditions */
			{
				double val = at(0, 0, 0) + 0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2]);
				double score = (Dxx + Dyy) * (Dxx + Dyy) / (Dxx * Dyy - Dxy * Dxy);
				double xn = x + b[0];
				double yn = y + b[1];
				double sn = s + b[2];

				vl_bool good =
					vl_abs_d(val) > tp &&
					score < (te + 1) * (te + 1) / te &&
					score >= 0 &&
					vl_abs_d(b[0]) < 1.5 &&
					vl_abs_d(b[1]) < 1.5 &&
					vl_abs_d(b[2]) < 1.5 &&
					xn >= 0 &&
					xn <= w - 1 &&
					yn >= 0 &&
					yn <= h - 1 &&
					sn >= s_min &&
					sn <= s_max;

				if (good)
				{
					k->o = f->o_cur;
					k->ix = x;
					k->iy = y;
					k->is = s;
					k->s = sn;
					k->x = xn * xper;
					k->y = yn * xper;
					k->sigma = f->sigma0 * pow(2.0, sn / f->S) * xper;
					++k;
				}

			} /* done checking */
		}	 /* next keypoint to refine */

		/* update keypoint count */
		f->nkeys = (int)(k - f->keys);
	}
}

// Normalizes descriptor using L1 Root Norm: L1-normalization followed by 
// an element-wise square root. This normalization is usually better than 
// standard L2-normalization.  See "Three things everyone should know to 
// improve object retrieval", Relja Arandjelovic and Andrew Zisserman, CVPR 2012.

VL_INLINE void L1Root_normalize(vl_sift_pix *begin, vl_sift_pix *end)
{
	vl_sift_pix *iter;
	vl_sift_pix norm = 0.0;

	for (iter = begin; iter != end; ++iter)
		norm += fabsf(*iter);
	norm += VL_EPSILON_F;

	for (iter = begin; iter != end; ++iter)
	{
		*iter /= norm;
		*iter = sqrtf(*iter);
	}
}
