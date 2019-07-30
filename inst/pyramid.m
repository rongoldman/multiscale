## -*- texinfo -*-
## @deftypefn{Function File} {p} = pyramid (@var{f}, @var{filter})
## Create a cell array containing a Gaussian pyramid decomposition othe field @var{f}. The @var{i+1}-th cell of @var{p} contains subsampling of the filtered @var{i}-th cell with the 2D filter @var{filter}. If @var{filter} is unsepecified, a 5x5 normalized Hamming filter is used.
## @seealso{filter2,hamming}
## @end deftypefn
function [p]=pyramid(f,fltr)
  if (nargin==1)
    fltr=hamming(5)*hamming(5)'/sum(hamming(5))^2;
 %   fltr=[1 1 1;1 1 1;1 1 1]/9;
  endif
  p=cell(1,max(round(log(size(f))/log(2))));
  p{1}=f;
  for i=2:length(p)
    p{i}=(filter2(fltr,p{i-1}(:,:),'same')./filter2(fltr,ones(size(p{i-1})),'same'))(1:2:end,1:2:end);
  endfor
endfunction
