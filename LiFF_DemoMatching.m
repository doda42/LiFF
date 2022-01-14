% LiFF_DemoMatching - Simple demonstration matching LiFF Light Field Features
% 
% This builds on LiFF_DemoSimple by locating features in two light fields, then matching them using
% vl_sift's vl_ubcmatch. Don't forget to run LiFF_PathSetup so that vl_feat's functions can be
% found by matlab.
% 
% See also LiFF_DemoSimple, LiFF_DemoColmapOut.m, LiFF_DemoFocalStack.m

% Part of LiFF Light Field Feature Toolbox v0.0.2
% Copyright (c) 2019-2022 Donald G. Dansereau

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

%---Rotate to make matching a bit more interesting---
LFSize = size(LF);
LF2 = zeros(LFSize, 'like', LF);
for( TIdx = 1:LFSize(1) )
	for( SIdx = 1:LFSize(1) )
		CurSlice = squeeze(LF(TIdx,SIdx,:,:));
		CurSlice = imrotate(CurSlice, 10, 'crop'); % not a good transform, LF2 is now misaligned
		LF2(TIdx,SIdx, :,:) = CurSlice;
	end
end

%---Find features and descriptors---
fprintf('Extracting features...\n');
tic
[f1,d1] = LiFF_ExtractFeatures( LF );
[f2,d2] = LiFF_ExtractFeatures( LF2 );
[matches, scores] = vl_ubcmatch(d1, d2);
toc

%---Display---
fprintf('Displaying matches...\n');
figure(1);
clf
Thumb1 = squeeze(LF(ceil(end/2),ceil(end/2),:,:));
Thumb2 = squeeze(LF2(ceil(end/2),ceil(end/2),:,:));
Thumb = cat(2,Thumb1,Thumb2);
imshow(Thumb);

NumMatches = size(matches,2); 
which=ceil(linspace(1,NumMatches,20)); % show 20 matches
x1 = f1(1,matches(1,which));
x2 = f2(1,matches(2,which)) + size(Thumb1,2);
y1 = f1(2,matches(1,which));
y2 = f2(2,matches(2,which));
 
hold on;
plot([x1; x2], [y1; y2]);

title(sprintf('LiFF -- %d matches (showing 20)', size(matches,2)))

fprintf('Done.\n');
