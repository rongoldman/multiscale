## -*- texinfo -*-
## @deftypefn{Function File} {f} = expand (@var{p}, @var{level}, @var{n}, @var{filter})
## Expand @var{n} times the @var{level} cell of pyramid @var{p} using the 2D filter @var{filter}. If @var{filter} is unsepecified, a 5x5 normalized Hamming filter is used.
## @seealso{pyramid,filter2,hamming}
## @end deftypefn
function [pout]=expand(p,j,n,fltr)
  if (!(j<=length(p) && j-n>=1)) 
    error('Index out of domain');
  endif
  if (nargin<4)
    fltr=hamming(5)*hamming(5)'/sum(hamming(5))^2;
  endif
  pout=p{j};
  for i=1:n
    pout=_expand1(p{j-i},pout,fltr);
  endfor
endfunction
