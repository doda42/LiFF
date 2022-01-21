% CIRCLE - Draws a circle.
%
% Usage: circle(c, r, n, col)
%
% Arguments:  c -  A 2-vector [x y] specifying the centre.
%             r -  The radius.
%             n -  Optional number of sides in the polygonal approximation.
%                  (defualt is 16 sides)
%           col -  optional colour, defaults to blue.

% Copyright (c) 1996-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Minor modifications for inclusion in LiFF Toolbox
% Copyright (c) 2019 Donald G. Dansereau

function circle(c, r, nsides, col, varargin)
nsides = LiFF_DefaultVal('nsides', 16);
col = LiFF_DefaultVal('col', [0,0,1]);

nsides = max(round(nsides),3); % Make sure it is an integer >= 3

a = [0:2*pi/nsides:2*pi];
line(r*cos(a)+c(1), r*sin(a)+c(2), 'color', col, varargin{:});
