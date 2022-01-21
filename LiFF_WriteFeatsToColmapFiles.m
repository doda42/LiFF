% LiFF_WriteFeatsToColmapFiles - Output features to files that can be loaded by COLMAP
% 
% Usage: 
%   LiFF_WriteFeatsToColmapFiles( OutFullPath, f,d )
%
% See LiFF_DemoColmapOut.m for example usage.
% 
% Inputs:
%   OutFullPath: path to output file
% 
%   f,d: features and feature descriptors to save; use LiFF_ConvertLiFFToColmap or
%   LiFF_ConvertVlToColmap to convert features and descriptors prior to calling
%   LiFF_WriteFeatsToColmapFiles.
% 
% See also LiFF_ReadFeatsFromColmapFiles.m, LiFF_ConvertLiFFToColmap.m, LiFF_ConvertVlToColmap.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function LiFF_WriteFeatsToColmapFiles( OutFullPath, f,d )

NumFeats = size(f,2);

OutFile = fopen(OutFullPath, 'wt');
fwrite( OutFile, sprintf('%d 128\n', NumFeats) );
for( iFeat = 1:NumFeats )
	CurFeat = f(:,iFeat);
	CurDesc = d(:,iFeat);
	CurFeatStr = sprintf('%g ', CurFeat);
	CurDescStr = sprintf('%d ', CurDesc);
	CurOutStr = sprintf('%s %s \n', CurFeatStr, CurDescStr);
	fwrite( OutFile, CurOutStr );
end
fclose( OutFile );
