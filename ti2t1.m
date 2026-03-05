function t1 = ti2t1(ti,tr)
%
% return t1-to-null from ti at specified tr
%

% arg checks
if ndims(ti)~=2 || ndims(tr)~=2
    error('ti and tr must be scalar, vector or matrix');
end

% broadcasting rules
t1 = zeros(size(ti.*tr));
ti = repmat(ti,size(t1)./size(ti));
tr = repmat(tr,size(t1)./size(tr));

% find the t1 with the given nulltime at the given tr
f = @(t1,ti,tr) t1.*(log(2)-log(1+exp(-tr./t1))) - ti;

% loop over dimensions
for k = 1:numel(t1)

    t1(k) = fzero(@(x)f(x,ti(k),tr(k)),ti(k)/log(2));

end
