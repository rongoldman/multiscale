## -*- texinfo -*-
## @deftypefn{Function File} {l} = lpyramid (@var{p}, @var{m}, @var{filter})
## Create a cell array containing a Laplacian pyramid based on the Gaussian pyramid @var{p}. Each cell of @var{l} contains the differance between a level of @var{p} and the next level of @var{p}, expanded and filtered with the 2D filter @var{filter}. If @var{filter} is unsepecified, a 5x5 normalized Hamming filter is used.
## @seealso{pyramid,filter2,hamming}
## @end deftypefn
function [l]=lpyramid(p,m,fltr)
  if(nargin<=3)
    #fltr=[0.25,0.5,0.25;0.5,1,0.5;0.25,0.5,0.25];
    fltr=hamming(5)*hamming(5)'/sum(hamming(5))^2;
  endif
  for i=length(p)-1:-1:1
    l{i}=p{i}-_expand1(p{i},p{i+1},fltr);
  endfor
endfunction
