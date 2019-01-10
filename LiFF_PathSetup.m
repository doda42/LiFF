% LiFF_PathSetup.m - add the LiFF light field feature toolbox to matlab's path
%
% It may be convenient to add this to your startup.m file.

% Part of LiFF Light Field Feature Toolbox v0.0.1
% Copyright (c) 2019 Donald G. Dansereau

function LiFF_PathSetup

if (~isdeployed)
	% Find the path to this script, and use it as the base path
	LiFFPath = fileparts(mfilename('fullpath'));
	
	fprintf('Adding path for LiFF Toolbox ');
	addpath( fullfile(LiFFPath) );
	
	fprintf('%s, done.\n', LiFF_ToolboxVersion);
end
