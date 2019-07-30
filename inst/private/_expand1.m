function p01=_expand1(p0,p1,fltr)
  p01=zeros(2*size(p1));
  if (min(size(p0(1:2:end,1:2:end))-size(p1))<0)
    size(p0(1:2:end,1:2:end)),size(p1)
  endif
  p01(1:2:end,1:2:end)=p1(:,:);
  p01(2:2:end,1:2:end)=p1(:,:);
  p01(1:2:end,2:2:end)=p1(:,:);
  p01(2:2:end,2:2:end)=p1(:,:);
  p01=resize(p01,size(p0));
  p01=filter2(fltr,p01,'same')./filter2(fltr,ones(size(p01)),'same');
endfunction
