% LiFF_ConvertLiFFToColmap - Prepares LiFF features for use by COLMAP
% 
% Usage: 
%   [f,d] = LiFF_ConvertLiFFToColmap( f,d )
% 
% This function prepares LiFF features for use by COLMAP, including correction to
% zero-based and 1/2-pixel-offset pixel indices and a change in ordering of the descriptors to UBC
% standard ordering.
% 
% For features extracted using VL_SIFT, use LiFF_ConvertVlToColmap.
%
% Inputs:
%   f: A list of features as detected by VL_SIFT
%   d: A matrix of L2-normalized floating point descriptors, as returned by VL_SIFT with the
%   'FloatDescriptors' option.
% 
% Outputs:
%    f,d: the converted features and descriptors ready for exporting to and use by COLMAP
% 
% See also:  LiFF_ConvertVlToColmap, LiFF_WriteFeatsToColmapFiles

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function [f,d] = LiFF_ConvertLiFFToColmap( f,d )

f(1:2,:) = f(1:2,:)-0.5;         % 0-based indexing, and pixel centers are at 0.5-pixel offsets
d = LiFF_ConvertToUBCDesc( d );  % UBC standard ordering is different from vl_feat standard

