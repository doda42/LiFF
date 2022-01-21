% LiFF_ConvertVlToColmap - Prepares VL_SIFT features for use by COLMAP
% 
% Usage: 
%   [f,d] = LiFF_ConvertVlToColmap( f,d )
% 
% This function prepares features extracted using VL_SIFT for use by COLMAP, including correction to
% zero-based and 1/2-pixel-offset pixel indices, a change in ordering of the descriptors to UBC
% standard ordering, and conversion of L2 normalization to L1 root normalization.
% 
% For LiFF features use LiFF_ConvertLiFFToColmap.
%
% The input includes a floating point descriptor, which can be obtained from VL_SIFT using the
% 'FloatDescriptors' parameter, as in 
%   [f,d] = vl_sift( Img, 'FloatDescriptors');
% 
% Inputs:
%   f: A list of features as detected by VL_SIFT
%   d: A matrix of L2-normalized floating point descriptors, as returned by VL_SIFT with the
%   'FloatDescriptors' option.
% 
% Outputs:
%    f,d: the converted features and descriptors ready for exporting to and use by COLMAP; the
%         descriptor will be in uint8 format as part of the conversion to L1 root norm.
% 
% See also:  LiFF_ConvertLiFFToColmap, LiFF_WriteFeatsToColmapFiles

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function [f,d] = LiFF_ConvertVlToColmap( f,d )

if( ~isfloat(d) )
	warning('Floating point descriptors are recommended when converting normalization methods; pass ''FloatDescriptors'' to vl_sift');
end

f(1:2,:) = f(1:2,:)-0.5;               % 0-based indexing, and pixel centers are at 0.5-pixel offsets
d = LiFF_ConvertToUBCDesc( d );      % UBC standard ordering is different from vl_feat standard
d = LiFF_ConvertL2ToL1RootNorm( d );   % use L1 root norm
