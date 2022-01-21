% LiFF_PathSetup.m - add the LiFF light field feature toolbox to matlab's path
%
% It may be convenient to add this to your startup.m file.
%
% Inputs:
%   SkipVLSetup (optional) : set this to true to skip calling vl_setup from the vl_sift that ships
%                            with LiFF. Useful when using your own version of vl_sift.

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function LiFF_PathSetup( SkipVLSetup )

if (~isdeployed)
	% Find the path to this script, and use it as the base path
	LiFFPath = fileparts(mfilename('fullpath'));
	
	fprintf('Adding path for LiFF Toolbox ');
	addpath( fullfile(LiFFPath) );

	% Optinally Call the bundled vl_setup to allow us to run vl_sift functions
	if( ~exist( 'SkipVLSetup', 'var' ) || SkipVLSetup==false )
		CurPath = pwd;
		cd( fullfile(LiFFPath, 'lib','vlfeat-0.9.21','toolbox') );
		vl_setup
		cd( CurPath );
	end
	
	fprintf('%s, done.\n', LiFF_ToolboxVersion);
end
