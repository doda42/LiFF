% LiFF_DemoFocalStack - Simple demonstration of light field focal stack
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
LF = LiFF_ReadESLF(InFile);
LF = LF(2:end-2,2:end-2,:,:,:);  % trim edge pixels, force odd view count
LF = single(LF);        % convert to float
LF = LF ./ max(LF(:));  % normalize

%---Find Focal Stack---
FocSteps = floor(size(LF,1)/2)*2 + 1;  % force to an odd number of steps
SlopeVec = linspace(-1, 1, FocSteps);
for( iChan=1:3 )
	FocStack(:,:,:,iChan) = LiFF_FocalStack( LF(:,:,:,:,iChan), SlopeVec );
end

%---Display---
FocStack = FocStack.^(1/2); % gamma
figure(1);
clf
for( iSlice = [1:size(FocStack,1), size(FocStack,1)-1:-1:1] )
	Thumb = squeeze(FocStack(iSlice,:,:,:));
	imshow(Thumb);
	drawnow
	pause(1/20)
end

