function varargout = argparse(varargin)
%VL_ARGPARSE Parse list of parameter-value pairs.
%   OPTS = VL_ARGPARSE(OPTS, ARGS) updates the structure OPTS based on
%   the specified parameter-value pairs ARGS={PAR1, VAL1, ... PARN,
%   VALN}. If a parameter PAR cannot be matched to any of the fields
%   in OPTS, the function generates an error.
%
%   Parameters that have a struct value in OPTS are processed
%   recursively, updating the individual subfields.  This behaviour
%   can be suppressed by using VL_ARGPARSE(OPTS, ARGS, 'nonrecursive'),
%   in which case the struct value is copied directly (hence deleting any
%   existing subfield existing in OPTS). A direct copy occurrs also if the
%   struct value in OPTS is a structure with no fields. The nonrecursive
%   mode is especially useful when a processing time is a concern.
%
%   One or more of the (PAR, VAL) pairs in the argument list can be
%   replaced by a structure; in this case, the fields of the structure
%   are used as paramater names and the field values as parameter
%   values. This behaviour, while orthogonal to structure-valued parameters,
%   is also disabled in the 'nonrecursive' mode.
%
%   A shorthand notation for substructs is available: the pair
%   ('field.subfield', value) creates a substruct, i.e. ('field',
%   struct('subfield', value)). It also works recursively (sub-sub-fields).
%
%   [OPTS, ARGS] = VL_ARGPARSE(OPTS, ARGS) copies any parameter in
%   ARGS that does not match OPTS back to ARGS instead of producing an
%   error. Options specified as structures are passed back as a list
%   of (PAR, VAL) pairs.
%
%   A further option is VL_ARGPARSE(OPTS, ARGS, 'merge'), which merges any
%   parameter in ARGS that do not match OPTS into OPTS. This essentially
%   turns off errors for unknown options, merging them anyway, and is
%   useful if OPTS does not necessarily define all possible default values.
%   The second output argument will always be empty. Another view is that
%   VL_ARGPARSE(OPTS1, OPTS2, 'merge') merges two structs into one.
%
%   Example::
%     The function can be used to parse a list of arguments
%     passed to a MATLAB functions:
%
%        function myFunction(x,y,z,varargin)
%        opts.parameterName = defaultValue ;
%        opts = vl_argparse(opts, varargin)
%
%     If only a subset of the options should be parsed, for example
%     because the other options are interpreted by a subroutine, then
%     use the form
%
%        [opts, varargin] = vl_argparse(opts, varargin)
%
%     that copies back to VARARGIN any unknown parameter.
%
%   See also: VL_HELP().
[varargout{1:nargout}] = vl_argparse(varargin{:});
