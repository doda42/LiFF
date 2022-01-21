% LiFF_ExtractFeatures - Detect and describe LiFF light field features
% 
% Usage: 
%   [F, D] = LiFF_ExtractFeatures( LF )
%   [F, D, FocStack] = LiFF_ExtractFeatures( LF, <OptionString>, <OptionVal>, ...)
%
% Detects LiFF keypoints and extracts their descriptors. This function closely mirrors the syntax
% and functionality of VL_SIFT, the VL_FEAT implementation of the SIFT feature detector and
% descriptor. For a demonstration see LiFF_DemoSimple.m.
%
% A full description of the LiFF feature detector and descriptor are available here:
%
%   [1] D. G. Dansereau, B. Girod, and G. Wetzstein, “LiFF: Light field features in scale and depth,” 
%       in Computer Vision and Pattern Recognition (CVPR), 2019. 
%       Paper and supplemental information are at https://roboticimaging.org/Tools/LiFF/
% 
% Inputs:
%   LF : A gray-scale 4D light field as single-precision floats. Index order is [t,s,v,u] with s,t
%   corresponding to horizontal and vertical subaperture indices, and u,v corresponding to
%   horisontal and vertical pixel indices.
% 
%   <OptionString> and <OptionVal>: optional arguments as described below.
% 
% Outputs:
%    F : detected features; Each column of F is a feature frame and has the format
%        [X; Y; S; TH; L], where X,Y are the pixel coordinates of the center of the center of the
%        frame, S is the scale, TH is the orientation in radians, and L is the index of the focal
%        stack slice contained the feature
%    D:  Each column of D is the descriptor of the corresponding frame in F. A descriptor is a 
%        128-dimensional vector of class UINT8. By default normalized using L1 Root normalization,
%        this can be overriden with the "L2Norm" option.
%    FocStack: Returns the focal stack constructed during the feature detection process, see
%    LiFF_FocalStack and LiFF_DemoFocalStack.
%
% Options:
%   FirstSlope, default -1:
%   LastSlope, default 1:
%     These set the range of slopes (depths) over which the feature search is executed
% 
%  NumSlopes:
%     Sets the number of slopes (depths) at which to compute the 4D slope/scale search space; the
%     default automatically covers the range -1 to 1 over as many steps as there are subaperture
%     views, rounded to the next largest odd number
% 
%   L2Norm:
%     Force L2 normalization of the descriptors, rather than the default L1 Root normalization. For
%     more on L1 Root normalization See "Three things everyone should know to improve object
%     retrieval", Relja Arandjelovic and Andrew Zisserman, CVPR 2012.
% 
%   FloatDescriptors:
%     Return descriptors in floating-point rather than uint8.
%
%   Octaves, default to maximum possible:
%     Set the number of octaves in the 4D slope/scale search space.
%
%   Levels, default 3:
%     Set the number of levels per octave of the 4D slope/scale search space.
%
%   FirstOctave, default 0:
%     Set the index of the first octave of the 4D slope/scale search space.
%
%   PeakThresh, default 0:
%     Set the peak selection threshold.
%
%   EdgeThresh, default 10:
%     Set the non-edge selection threshold.
%
%   NormThresh, default -inf:
%     Set the minimum l2-norm of the descriptors before normalization. Descriptors below the
%     threshold are set to zero.
%
%   Magnif, default 3:
%     Set the descriptor magnification factor. The scale of the keypoint is multiplied by this
%     factor to obtain the width (in pixels) of the spatial bins. For instance, if there are there
%     are 4 spatial bins along each spatial direction, the ``side'' of the descriptor is
%     approximatively 4 * MAGNIF.
%
%   WindowSize, default 2:
%     Set the variance of the Gaussian window that determines the descriptor support. It is
%     expressend in units of spatial bins.
%
%   Frames:
%     If specified, set the frames to use (bypass the detector). If frames are not passed in order
%     of increasing scale, they are re-orderded.
%
%   Orientations:
%     If specified, compute the orientations of the frames overriding the orientation specified by
%     the 'Frames' option.
%
%   Verbose:
%     If specfified, be verbose. May be repeated to increase the verbosity level.
% 
% Example:
%      [f,d] = LiFF_ExtractFeatures( LF, 'FirstOctave', -1 );
% 
% See also:  LiFF_DemoSimple, LiFF_DemoFocalStack, LiFF_DemoColmapOut, LiFF_FocalStack

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau
