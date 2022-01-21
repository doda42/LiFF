% LiFF_RGB2Gray.m - Convert colour light field to grayscale
% 
% Usage: 
%   LF = LiFF_RGB2Gray( LF )
%
% See LiFF_DemoSimple.m for example usage.
%
% This function applies the MATLAB rgb2gray function across the light field, then applies a gamma
% correction factor of 1/2.
% 
% See also LiFF_DemoSimple.m, LiFF_DemoFocalStack.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function LF = LiFF_RGB2Gray( LF )

LFSize = size(LF);

% Apply matlab's built-in perceptually based rgb2gray
LF = reshape(LF, prod(LFSize([1,3])), prod(LFSize([2,4])), 3);
LF = rgb2gray(LF);
LF = reshape(LF, LFSize(1:4));

% Apply gamma 
% Standard gamma is 1/2.2, but 1/2 is much faster
LF = LF.^(1/2);
