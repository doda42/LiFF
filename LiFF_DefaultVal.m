% LiFF_DefaultVal - Convenience function to set up default parameter values
% 
% Usage: 
% 
%   Var = LiFF_DefaultVal( Var, DefaultVal )
% 
% 
% This provides an elegant way to establish default parameter values. 
%
% This function is adapted from LFDefaultVal from the Light Field Toolbox.
%
% Inputs:
% 
%   Var: string giving the name of the parameter
%   DefaultVal: default value for the parameter
%
% 
% Outputs:
% 
%   Var: if the parameter already existed, the output matches its original value, otherwise the
%        output takes on the specified default value
% 
% Example: 
% 
%   clearvars
%   ExistingVar = 42;
%   ExistingVar = LiFF_DefaultVal( 'ExistingVar', 3 )
%   OtherVar = LiFF_DefaultVal( 'OtherVar', 3 )
% 
%   Results in :
%       ExistingVar =
%           42
%       OtherVar =
%            3

% Part of LiFF Light Field Feature Toolbox
% Copyright (c) 2019 Donald G. Dansereau


function Var = LiFF_DefaultVal( Var, DefaultVal )

CheckIfExists = sprintf('exist(''%s'', ''var'') && ~isempty(%s)', Var, Var);
VarExists = evalin( 'caller', CheckIfExists );

if( ~VarExists )
    Var = DefaultVal;
else
    Var = evalin( 'caller', Var );
end

end