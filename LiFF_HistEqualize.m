% LiFF_HistEqualize - histogram-based contrast adjustment
%
% Usage: 
%
%     LF = LiFF_HistEqualize( LF )
%     LF = LiFF_HistEqualize( LF, [], LowerBound, UpperBound )
%     LF = LiFF_HistEqualize(LF, Cutoff_percent, [], [], Precision)
% 
% This function performs contrast adjustment on images based on a histogram of intensity. Colour
% images are handled by building a histogram of the `value' channel in the HSV colour space.
% Saturation points are computed based on a desired percentage of saturated white / black pixels.
%
% This function is adapted from the Light Field Toolbox function LFHistEqualize.
% 
% All parameters except LF are optional. Pass an empty array "[]" to omit a parameter.
%
% Inputs:
% 
%     LF : a light field to adjust. 
% 
%          Including a weight channel will result in zero-weight pixels being ignored. Weight should
%          be included as a fourth colour channel -- e.g. an image of size [Nl,Nk,4].
%
%  Optional inputs:
% 
%             Cutoff_percent : controls how much of the histogram gets saturated on both the high
%                              and lower ends. Higher values result in more saturated black and
%                              white pixels. The default value is 0.1%.
%
%     LowerBound, UpperBound : bypass the histogram computation and instead directly saturate the
%                              input image at the prescribed limits. If no bounds are provided, they
%                              are computed from the histogram.
% 
%                  Precision : 'single' or 'double'
%
% 
% Outputs:
% 
%                         LF : the contrast-adjusted image, in floating point values scaled between
%                              0 and 1.
%
%     LowerBound, UpperBound : the saturation points used, on the intensity scale of the input after
%                              converting to the requested floating point precision. These are
%                              useful in applying the same contrast adjustment to a number of light
%                              fields, by saving and passing these bounds on subsequent function
%                              calls.
%
% See LiFF_DemoColmapOut.m for example usage.

% Adapted from Light Field Toolbox 0.4
% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function [LF, LowerBound, UpperBound] = ...
    LFHistEqualize(LF, Cutoff_percent, LowerBound, UpperBound, Precision)

%---Defaults---
if( isfloat(LF) )
	DefaultPrecision = class(LF);
else
	DefaultPrecision = 'single';
end
Cutoff_percent = LiFF_DefaultVal( 'Cutoff_percent', 0.1 );
Precision = LiFF_DefaultVal( 'Precision', DefaultPrecision );
MaxBlkSize = LiFF_DefaultVal( 'MaxBlkSize', 10e6 ); % slice up to limit memory utilization

%---Make sure LF is floating point, as required by hist function---
if( ~strcmp(class(LF), Precision) )
    LF = cast(LF, Precision);
end

%---Flatten to a single list of N x NChans entries---
LFSize = size(LF);
NDims = numel(LFSize);
NChans = LFSize(end);
LF = reshape(LF, [prod(LFSize(1:NDims-1)), NChans]);

%---Strip out and use the weight channel if it's present---
if( NChans == 4 || NChans == 2 )
    % Weight channel found, strip it out and use it
    LFW = LF(:,end);
    LF = LF(:,1:end-1);
    ValidIdx = find(LFW > 0);
else
    LFW = [];
    ValidIdx = ':';
end

%---Convert colour images to HSV, and operate on value only---
if( size(LF,2) == 3 )
    LF = LF - min(LF(ValidIdx));
    LF = LF ./ max(LF(ValidIdx));
    LF = max(0,min(1,LF));
    
    for( UStart = 1:MaxBlkSize:size(LF,1) )
        UStop = UStart + MaxBlkSize - 1;
        UStop = min(UStop, size(LF,1));

        % matlab 2014b's rgb2hsv requires doubles
        LF(UStart:UStop,:) = cast(rgb2hsv( double(LF(UStart:UStop,:)) ), Precision);  
    end
    LF_hsv = LF(:,1:2);
    LF = LF(:,3); % operate on value channel only
else
    LF_hsv = [];
end

%---Compute bounds from histogram---
if( ~exist('LowerBound','var') || ~exist('UpperBound','var') || ...
        isempty(LowerBound) || isempty(UpperBound))
    [ha,xa] = hist(LF(ValidIdx), 2^16);
    ha = cumsum(ha);
    ha = ha ./ max(ha(:));
    Cutoff = Cutoff_percent / 100;
    LowerBound = xa(find(ha >= Cutoff, 1, 'first'));
    UpperBound = xa(find(ha <= 1-Cutoff, 1, 'last'));
end
if( isempty(UpperBound) )
    UpperBound = max(LF(ValidIdx));
end
if( isempty(LowerBound) )
    LowerBound = min(LF(ValidIdx));
end
clear ValidIdx

%---Apply bounds and rescale to 0-1 range---
LF = max(LowerBound, min(UpperBound, LF));
LF = (LF - LowerBound) ./ (UpperBound - LowerBound);

%---Rebuild hsv, then rgb colour light field---
if( ~isempty(LF_hsv) )
    LF = cat(2,LF_hsv,LF);
    clear LF_hsv
    
    for( UStart = 1:MaxBlkSize:size(LF,1) )
        UStop = UStart + MaxBlkSize - 1;
        UStop = min(UStop, size(LF,1));
        
        LF(UStart:UStop,:) = hsv2rgb(LF(UStart:UStop,:));
    end
end

%---Return the weight channel---
if( ~isempty(LFW) )
    LF(:,end+1) = LFW;
end

%---Unflatten---
LF = reshape(LF, [LFSize(1:NDims-1),NChans]);
