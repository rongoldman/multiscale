## Copyright (C) 2016-2017 Carlo de Falco
## Copyright (C) 2016 Francesco Faccio <francesco.faccio@mail.polimi.it>
## Copyright (C) 2013-2016 Roberto Porcu' <roberto.porcu@polimi.it>
## Copyright (C) 2006-2012 Thomas Treichl <treichl@users.sourceforge.net>
##
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {@var{mgstruct} =} mgset ()
## @deftypefnx {} {@var{mgstruct} =} mgset (@var{"field1"}, @var{value1}, @var{"field2"}, @var{value2}, @dots{})
## @deftypefnx {} {@var{mgstruct} =} mgset (@var{oldstruct}, @var{"field1"}, @var{value1}, @var{"field2"}, @var{value2}, @dots{})
## @deftypefnx {} {@var{mgstruct} =} mgset (@var{oldstruct}, @var{newstruct})
## @deftypefnx {} {} mgset ()
##
## Create or modify an Multigrid options structure.
##
## When called with no input argument and one output argument, return a new Multigrid
## options structure that contains all possible fields initialized to their
## default values.  If no output argument is requested, display a list of
## the common Multigrid solver options along with their default value.
##
## If called with name-value input argument pairs @var{"field1"},
## @var{"value1"}, @var{"field2"}, @var{"value2"}, @dots{} return a new
## Multigrid options structure with all the most common option fields
## initialized, @strong{and} set the values of the fields @var{"field1"},
## @var{"field2"}, @dots{} to the values @var{value1}, @var{value2},
## @dots{}.
##
## If called with an input structure @var{oldstruct} then overwrite the
## values of the options @var{"field1"}, @var{"field2"}, @dots{} with
## new values @var{value1}, @var{value2}, @dots{} and return the
## modified structure.
##
## When called with two input Multigrid options structures @var{oldstruct} and
## @var{newstruct} overwrite all values from the structure
## @var{oldstruct} with new values from the structure @var{newstruct}.
## Empty values in @var{newstruct} will not overwrite values in
## @var{oldstruct}.
##
## The most commonly used Multigrid options, which are always assigned a value
## by @qcode{mgset}, are the following:
##
## @table @asis
## @item Tol
## Error tolerance. If set to 0.0 the stopping criteria is not applied.
##
## @item FrontIterations
## Number of iterations of the front relaxation
##
## @item BackIterations
## Number of iterations of the back relaxation
##
## @item InnerIterations
## Number of iterations of solving the coarsest grid
##
## @item MaxIterations
## Maximum number of iterarions at the finest grid. If Tol is used and is 
## not reached, a failed status is reported.
##
## @item FrontRelaxation
## Relaxation scheme applied before coarsening. Accepted values include
## @qcode{"J"}, @qcode{"RB"}, @qcode{"Z1"}, @qcode{"Z2"}.
##
## @item BackRelaxation
## Relaxation scheme applied after refining. Accepted values include
## @qcode{"J"} for Jaccobi, @qcode{"RB"} for red-black Gauss-Sidel (default), 
## @qcode{"Z1"} for zebra on 1st dimension, @qcode{"Z2"} for zebra on 2nd dimension.
##
## @item NormControl
## Norm to control error relative to the  of the solution, 2-norm and
## absolute value (inf-norm) are implemented.
##
## @item Verbose
## Should the program provide output
##
## @end table
##
## Field names that are not in the above list are also accepted and
## added to the result structure.
##
## @seealso{mgget}
## @end deftypefn

function mgstruct = mgset (varargin)

  persistent p;

  if (isempty (p))
    p = inputParser ();
    p.addParameter ("Tol", []);
    p.addParameter ("FrontIterations", []);
    p.addParameter ("BackIterations", []);
    p.addParameter ("InnerIterations", []);
    p.addParameter ("MaxIterations", []);
    p.addParameter ("FrontRelaxation", []);
    p.addParameter ("BackRelaxation", []);
    p.addParameter ("NormControl", []);
    p.addParameter ("Verbose", []);
    p.KeepUnmatched = true;
  endif

  if (nargin == 0 && nargout == 0)
    print_options ();
  else
    p.parse (varargin{:});
    mgstruct = p.Results;
    mgstruct_extra = p.Unmatched;

    xtra_fields = fieldnames (mgstruct_extra);
    if (! isempty (xtra_fields))
      ## Merge extra fields into existing mgstruct
      for fldname = sort (xtra_fields.')
        fldname = fldname{1};
        warning ("Octave:invalid-input-arg",
                 "mgset: unknown option \"%s\"\n", fldname);
        mgstruct.(fldname) = mgstruct_extra.(fldname);
      endfor
    endif

  endif

endfunction

## function to print all possible options
function print_options ()

  disp ("List of the most common Multigrid solver options.");
  disp ("Default values are in square brackets.");
  disp ("");
  disp ('             Tol:  non-negative scalar, [0.0]');
  disp (' FrontIterations:  positive integer,      [1]');
  disp ('  BackIterations:  positive integer,      [1]');
  disp (' InnerIterations:  positive integer,      [1]');
  disp ('   MaxIterations:  positive integer,     [20]');
  disp (' FrontRelaxation:  switch, {"J",["RB"], "Z1", "Z2"}');
  disp ('  BackRelaxation:  switch, {"J",["RB"], "Z1", "Z2"}');
  disp ('     NormControl:  switch, {[2], inf}');
  disp ('         Verbose:  switch, {true,[false]}');

endfunction


%!demo
%! ## A new Multigrid options structure with default values is created.
%!
%! mgoptA = mgset ();

%!demo
%! ## A new Multigrid options structure with manually set options
%! ## for "Tol" is created.
%!
%! mgoptB = mgset ("Tol", 1e-2);

%!demo
%! ## A new Multigrid options structure is created from mgoptB with
%! ## a modified value for option "NormControl".
%!
%! mgoptB = mgset ("Tol", 1e-2);
%! mgoptC = mgset (mgoptB, "NormControl", inf);

%!test
%! mgoptA = mgset ();
%! assert (isstruct (mgoptA));
%! assert (numfields (mgoptA), 8);
%! assert (all (structfun ("isempty", mgoptA)));

%!shared mgoptB, mgoptC
%!test
%! mgoptB = mgset ("TOL", 1e-2, "maxiterations", 6);
%! assert (mgoptB.Tol, 1e-2);  # Check canonicalization of name
%! assert (mgoptB.MaxIterations, 6);

%!test
%! mgoptC = mgset (mgoptB, "NormControl", 2);
%! assert (mgoptC.Tol, 1e-2);       # check values from first struct copied
%! assert (mgoptC.NormControl, 2);  # check new values override old ones

%!test
%! mgoptD = mgset (mgoptB, mgoptC);
%! assert (mgoptD, mgoptC);

## Test custom user-defined option
%!test
%! warning ("off", "Octave:invalid-input-arg", "local");
%! mgopt = mgset ("NewtonTol", 3);
%! assert (mgopt.NewtonTol, 3);

## FIXME: Add an inexact match option once it is available in inputParser.
## See bug #49364.
## %!warning <no exact match for 'max'.  Assuming 'MaxIterations'> mgset ("max", 1);

## Test input validation
%!error <argument 'OPT1' is not a valid parameter> mgset ("opt1")
%!error mgset (1, 1)
%!error <argument 'OPT1' is not a valid parameter> mgset (mgset (), "opt1")
%!error mgset (mgset (), 1, 1)
