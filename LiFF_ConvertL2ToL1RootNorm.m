% LiFF_ConvertL2ToL1RootNorm - Converts L2-normalized descriptors to L1-root-normalized
% 
% Usage: 
%   d = LiFF_ConvertL2ToL1RootNorm( d )
% 
% This is useful when working with other feature detectors, e.g. VL_SIFT, that do not directly
% support L1 root normalization.  This is not necessary for LiFF toolbox features, which are L1 root
% normalized by default.  This function is used by LiFF_ConvertVlToColmap to prepare VL_SIFT
% features for use by COLMAP.  
%
% The input is a floating point descriptor, which can be obtained from VL_SIFT using the
% 'FloatDescriptors' parameter, as in 
%   [f,d] = vl_sift( Img, 'FloatDescriptors');
% 
% Inputs:
%   d: A matrix of L2-normalized floating point descriptors, as returned by vl_sift with the
%   'FloatDescriptors' option.
% 
% Outputs:
%    d: a matrix of descriptors renormalized using the L1 root norm and converted to uint8 format
% 
% See also:  LiFF_ConvertVlToColmap

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function d = LiFF_ConvertL2ToL1RootNorm( d )

if( ~isfloat(d) )
	warning('Floating point descriptors are recommended when converting normalization methods');
end

d = sqrt(d); 
k = 512 * sqrt(1./sum(d.^2));
d = d .* k;

d = uint8(floor(d)); 
