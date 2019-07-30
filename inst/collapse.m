## -*- texinfo -*-
## @deftypefn{Function File} {f} = collapse (@var{lp})
## @deftypefnx{Function File} {f} = collapse (@var{fcn}, @var{lp1}, @var{lp2}, @dots{})
## @deftypefnx{Function File} {f} = collapse (@dots{}, @var{filter})
## Collapse a Laplacian pyramid @var{lp} by expanfing each level to the first and adding all the expanded levels. If the arguments are a function handle @var{fcn} followed by one or more pyramids of the same size, then each expanded level is passed as an argument to @var{fcn} prior to addition of all levels. If @var{filter} is unsepecified, a 5x5 normalized Hamming filter is used.
## @seealso{pyramid, lpyramid, expand, filter2, hamming}
## @end deftypefn
function [f]=collapse(varargin)
for i=1:nargin
 typeinfo(varargin{i})
endfor
if (isa(varargin{1},'cell'))
 if (nargin==2 && isa(varargin{2},'matrix'))
  fltr=varargin{2};
 else
  fltr=hamming(5)*hamming(5)'/sum(hamming(5))^2;
 endif
 lp=varargin{1};
 f=zeros(size(lp{1}));
 for i=1:length(lp)
   f+=expand(lp,i,i-1,fltr);
 endfor
elseif (isa(varargin{1},'function_handle'))
 fnc=varargin{1};
 if (isa(varargin{end},'matrix'))
  fltr=varargin{end};
 else
  fltr=hamming(5)*hamming(5)'/sum(hamming(5))^2;
 endif
 f=zeros(size(varargin{2}{1}));
 for j=1:length(varargin{2})
  tmp=cell(1,nargin-1);
  for i=2:nargin
   tmp{i-1}=expand(varargin{i},j,j-1,fltr);
  endfor
  f+=fnc(tmp);
 endfor
endif

endfunction
