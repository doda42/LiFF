/*
FocalStack.c

C implementation of light field focal stack
See LiFF_FocalStack.m

Part of LiFF Light Field Feature Toolbox
Copyright (c) 2019 Donald G. Dansereau
*/


#include <mex.h>
#include <math.h>
#include <assert.h>

//---
void DoFocalStack( const float* pLF, const mwSize* pLFSize, const double* pSlopeVec, int NumSlopes, float* pOutStack, const mwSize* pStackSize )
{
	float TCent = (pLFSize[0]-1)/2.0f;
	float SCent = (pLFSize[1]-1)/2.0f;

	for( int iSlope=0; iSlope<pStackSize[0]; ++iSlope )
	{
		float CurSlope = pSlopeVec[iSlope];
		for( int VIdx=0; VIdx<pLFSize[2]; ++VIdx )
		{
			for( int UIdx=0; UIdx<pLFSize[3]; ++UIdx )
			{
				float CurVal = 0;
				int CurCount = 0;
				
				float CurV = VIdx + -TCent * CurSlope;
				for( int TIdx=0; TIdx<pLFSize[0]; ++TIdx )
				{
					int CurVIdx = (int)(CurV+0.5f);
					CurV += CurSlope;
					if( CurVIdx < 0 || CurVIdx >= pLFSize[2] )
					{
						continue;
					}

					float CurU = UIdx + -SCent * CurSlope;
					for( int SIdx=0; SIdx<pLFSize[1]; ++SIdx )
					{
						int CurUIdx = (int)(CurU+0.5f);
						CurU += CurSlope;
						if( CurUIdx < 0 || CurUIdx >= pLFSize[3] )
						{
							continue;
						}
						
						CurVal = CurVal + pLF[TIdx + SIdx*pLFSize[0] + CurVIdx*pLFSize[0]*pLFSize[1] + CurUIdx*pLFSize[0]*pLFSize[1]*pLFSize[2] ];
						++CurCount;
					}
				}
				CurVal /= CurCount; // normalize to reduce edge effects
				pOutStack[iSlope + VIdx*pStackSize[0] + UIdx*pStackSize[0]*pStackSize[1]] = CurVal;
			}
		}
	}
}

