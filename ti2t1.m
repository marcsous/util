function t1 = ti2t1(ti,tr)
%
% return t1 to null from ti at specified tr
%

if numel(tr)==1 && numel(ti)>1
    tr = repmat(tr,size(ti));
elseif numel(tr)>1 && numel(ti)==1
    ti = repmat(ti,size(tr));
elseif numel(tr)~=numel(ti)
    error('ti and tr incompatible');
end

for k = 1:numel(ti)

    f = @(t1)t1.*(log(2)-log(1+exp(-tr(k)./t1))) - ti(k);

    t1(k) = fzero(f,ti(k)/log(2));

end

t1 = reshape(t1,size(ti));