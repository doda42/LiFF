% LiFF_ConvertToUBCDesc - Reorder a LiFF or VL_SIFT descriptor to UBC standard ordering
% 
% Usage: 
%   [d] = LiFF_ConvertToUBCDesc( d )
% 
% LiFF and VL_SIFT use the same ordering for descriptors. This function converts this ordering to
% the UBC standard ordering used internally by COLMAP. This is helpful for making direct comparisons
% with COLMAP-generated features.
% 
% Inputs:
%   d: A matrix of descriptors in LiFF/VL_SIFT standard ordering
% 
% Outputs:
%   d: the re-ordered descriptors
% 
% See also:  LiFF_ConvertLiFFToColmap, LiFF_WriteFeatsToColmapFiles

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

% todo[optimization]: vectorizing this code should make it faster

function UbcDescriptors = LiFF_ConvertToUBCDesc( VlFeatDescriptors )
q = [0, 7, 6, 5, 4, 3, 2, 1];

UbcDescriptors = zeros(size(VlFeatDescriptors), 'like', VlFeatDescriptors);

for( n = 1:size(VlFeatDescriptors,2) )
	for( i = 0:4-1 )
		for( j = 0:4-1 )
			for( k = 0:8-1 )
				NewVal = VlFeatDescriptors(8 * (j + 4 * i) + k +1, n);
				UbcDescriptors(8 * (j + 4 * i) + q(k+1) +1, n) = NewVal;
			end
		end
	end
end