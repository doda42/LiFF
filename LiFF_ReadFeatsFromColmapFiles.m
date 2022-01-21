% LiFF_ReadFeatsFromColmapFiles - Read features from COLMAP files
% 
% Usage: 
%   [f,d] = LiFF_ReadFeatsFromColmapFiles( InFullPath )
% 
% Inputs:
%   OutFullPath: path to output file
% 
% Outputs:
%   f,d: features and feature descriptors loaded from file; feature locations are automatically
%   converted to 1-based indexing without 1/2 pixel offsets; descriptor ordering is not modified.
% 
% See also LiFF_WriteFeatsToColmapFiles.m, LiFF_ConvertLiFFToColmap.m, LiFF_ConvertVlToColmap.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau


function [f,d] = LiFF_ReadFeatsFromColmapFiles( InFullPath )

InFile = fopen(InFullPath, 'rt');
CurLine = fgetl(InFile);
NumFeats = sscanf(CurLine, '%d', 1);

for( iFeat = 1:NumFeats )
	CurLine = fgetl(InFile);
	CurFeat = sscanf(CurLine, '%g');
	CurDesc = CurFeat(6:end);
	CurFeat = CurFeat(1:5);

	f(:,iFeat) = CurFeat;
	d(:,iFeat) = CurDesc;
end
fclose( InFile );

f(1:2,:) = f(1:2,:)+0.5; % 0-based indexing, pixel centers are at 0.5 offsets

d = uint8(d);
