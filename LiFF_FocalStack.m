% LiFF_FocalStack - Compute focal stack from light field
% 
% Usage: 
%   FocStack = LiFF_FocalStack( LF, SlopeVec )
%
% Computes a focal stack, one slice for each of the entries in SlopeVec. See LiFF_DemoFocalStack.
% 
% The focal stack is computed by the "shift-and-add" method, in which subaperture images are shifted
% according to a given slope, then added together.  Interpolation is nearest-neighbour, and
% upsampling is not implemented. Values are normalized to reduce edge effects.
% 
% Inputs:
%   LF : A gray-scale 4D light field as single-precision floats. Index order is [t,s,v,u] with s,t
%   corresponding to horizontal and vertical subaperture indices, and u,v corresponding to
%   horisontal and vertical pixel indices.
% 
%   SlopeVec: A list of slopes at which to compute each focal stack slice.
% 
% Outputs:
%    FocStack: A focal stack, index order [slope,v,u].
% 
% Example:
%    FocStack = LiFF_FocalStack( LF, linspace(-1,1, 15) );
%    imshow( squeeze(FocStack(1,:,:)) );
% 
% See also:  LiFF_DemoFocalStack

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

% todo: add more interpolation options; compare vs fast focal stack algorithms; add upsampling
% options;