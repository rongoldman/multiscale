## -*- texinfo -*-
## @deftypefn{Function File} {} mid (@var{X}, @var{DIM})
## If @var{X} is a vector of length N, @code{diff (X)} is the vector of averages:
## @example 
## @group 
## 0.5(@var{X}(2) + @var{X}(1)), @dots{}, 0.5(@var{X}(n) + @var{X}(n-1)) 
## @end group 
## @end example.
## If @var{X} is a matrix, @code{mid (X)} is the matrix of column averages along the first non-singleton dimension.
##
## The dimension along which to take the averages can be explicitly stated with the optional variable DIM. 
##
## @end deftypefn

function y=mid(varargin)
  dims=size(varargin{1});
  if (length(varargin)>1) 
    dim=varargin{2};
  else
    dim=find(dims>1)(1);
  endif
  dims(dim)-=1;
  y=resize(shift(varargin{1},-1,dim)+varargin{1},dims)/2;
endfunction
