% LiFF_DefaultVal - Convenience function to set up default parameter values
% 
% Usage: 
% 
%   Var = LFDefaultVal( Var, DefaultVal )
% 
% 
% This provides an elegant way to establish default parameter values. See LFDefaultField for setting
% up structs with default field values.
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
%   ExistingVar = LFDefaultVal( 'ExistingVar', 3 )
%   OtherVar = LFDefaultVal( 'OtherVar', 3 )
% 
%   Results in :
%       ExistingVar =
%           42
%       OtherVar =
%            3

% Part of LiFF Light Field Feature Toolbox v0.0.1
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