% make_win64 - Build mex binaries for Windows 64
% 
% Usage: 
%   make_win64
% 
% See also make_linux.m

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau

mex -I'./lib/vlfeat-0.9.21' COPTIMFLAGS='-O3 -DNDEBUG' LiFF_FocalStack.c FocalStack.c
mex -I'./lib/vlfeat-0.9.21' COPTIMFLAGS='-O3 -DNDEBUG' -L'./lib/vlfeat-0.9.21/bin/win64' -l'vl' LiFF_ExtractFeatures.c FocalStack.c
