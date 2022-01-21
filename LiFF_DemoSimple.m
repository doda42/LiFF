% LiFF_DemoSimple - Simple demonstration of LiFF Light Field Features
% 
% The demo loads a light field, converts it to grayscale, locates features, and visualizes them
% using colour to indicate the slope at which each feature is identified. 
%
% Change the input file in the "Tweakables" section at the top of the script.
% 
% See also LiFF_DemoColmapOut.m, LiFF_DemoFocalStack.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

clearvars

%---Tweakables---
InFile = 'SampleScenes/Plant.eslf.jpg';

%---Load---
fprintf('Loading light field and converting to grayscale...\n');
LF = LiFF_ReadESLF(InFile);
LF = LF(2:end-2,2:end-2,:,:,:); % remove pixels on lenslet borders, 2 on the bot/right, 1 on top/left
LF = single(LF);        % convert to float
LF = LF ./ max(LF(:));  % normalize
LF = LiFF_RGB2Gray(LF); % convert to grayscale

%---Find features and descriptors---
fprintf('Extracting features...\n');
tic
[f,d] = LiFF_ExtractFeatures( LF );
toc

%---Display---
fprintf('Displaying features...\n');
figure(1);
clf
Thumb = squeeze(LF(ceil(end/2),ceil(end/2),:,:));
imshow(Thumb);

MinSlope = min(f(5,:));
SlopeRange = max(1, max(f(5,:)) - MinSlope);
for( i=1:size(f,2) )
	CurFeat = f(:, i);
	Color = (CurFeat(5)-MinSlope)/SlopeRange;
	RGBColor = [1-Color, 1-2*abs(Color-0.5), Color];
	circle( [CurFeat(1), CurFeat(2)], CurFeat(3), [],RGBColor, 'linewidth', 2 );
end
title(sprintf('LiFF -- %d features', size(f,2)))

fprintf('Done.\n');
