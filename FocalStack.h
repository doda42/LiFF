/*
FocalStack.h

Interface for C implementation of light field focal stack
See LiFF_FocalStack.m

Part of LiFF Light Field Feature Toolbox
Copyright (c) 2019 Donald G. Dansereau
*/

#ifndef FOCALSTACK_H
#define FOCALSTACK_H

void DoFocalStack( const float* pLF, const mwSize* pLFSize, const double* pSlopeVec, int NumSlopes, float* pOutStack, const mwSize* pStackSize );

#endif
