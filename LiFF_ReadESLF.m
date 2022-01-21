% LiFF_ReadESLF.m - Read ESLF-encoded light fields
% 
% Usage: 
%   LF = LiFF_ReadESLF( FName )
%   LF = LiFF_ReadESLF( FName, LensletSize_pix, HasAlpha )
%
% See LiFF_DemoSimple.m for example usage.
%
% This function reads the ESLF files commonly used to store light fields, e.g. those at
% lightfields.stanford.edu.
% 
% Inputs:
%   FName: path to input file
% 
%   LensletSize_pix, default 14: number of pixels per lenslet
% 
%   HasAlpha, default false: if true, a weight channel is added to the loaded light field indicating
%   which subaperture views are empty
% 
% Outputs:
%   LF: a 4D light field with index order [t,s,v,u,c], where s,t are horizontal and vertical
%   subaperture index; u,v are horizontal and vertical pixel index; and c is colour channel.
%
% See also LiFF_DemoSimple.m, LiFF_DemoFocalStack.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function LF = LiFF_ReadESLF( FName, LensletSize_pix, HasAlpha )

LensletSize_pix = LiFF_DefaultVal('LensletSize_pix',[14,14]);
HasAlpha = LiFF_DefaultVal('HasAlpha',false);
NChans = 3;

if( HasAlpha )
	[Img,ColorMap, Alpha] = imread( FName ); % note: requesting Alpha sets bg to black
	Img(:,:,4) = Alpha;
	NChans = NChans + 1;
else
	[Img,ColorMap] = imread( FName );
end

if( ~isempty(ColorMap) )  % indexed colour image
	Img = ind2rgb(Img, ColorMap); % todo[optimization] may wish to convert to uint8
end

ImgSize = size(Img(:,:,1));

LFSize(3:4) = ceil(ImgSize./LensletSize_pix);
PadSize = LFSize(3:4) .* LensletSize_pix;
PadAmt = PadSize-ImgSize;

LFSize(1:2) = PadSize ./ LFSize(3:4);
LFSize(5) = NChans;

ImgPad = padarray(Img, PadAmt, 0,'post');

LF = reshape(ImgPad, LFSize([1,3,2,4,5]));
LF = permute(LF, [1,3,2,4,5]);
