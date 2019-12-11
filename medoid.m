function y = medoid(x,dim,nanflag)
%y = medoid(x,dim,nanflag)
%
% Calculates medoid along a specified dimension.
%
%    MIN { SUM_i |y - x_i| } for y in x
%
% The medoid is similar to the median but always
% returns an element of x. It's also well defined
% for complex arguments, unlike the median.
%
%%

szx = size(x);

% error handling
if any(szx==0)
    error(' medoid: x cannot be an empty matrix');
end
if nargin<2
    [~,dim] = max(szx>1);
else
    validateattributes(dim,{'numeric'},{'positive','scalar','integer'},'','dim');
    if dim>numel(szx); y = x; return; end
end
if nargin<3 || isequal(nanflag,'includenan')
    omitnan = false;
elseif isequal(nanflag,'omitnan')
    omitnan = true;
else
    error('nanflag option not recognized.');
end

% size of return argument
szy = szx; szy(dim) = 1;

% evaluate sum(|d|) for all y (brute force)
if isa(x,'gpuArray')
    sumd = zeros(szx,'gpuArray');
else
    sumd = zeros(szx);
end

for k = 2:szx(dim)

    % form abs difference
    switch dim
        case 1; y = x(k-1,:,:,:,:); z = x(k:end,:,:,:,:);
        case 2; y = x(:,k-1,:,:,:); z = x(:,k:end,:,:,:);
        case 3; y = x(:,:,k-1,:,:); z = x(:,:,k:end,:,:);
        case 4; y = x(:,:,:,k-1,:); z = x(:,:,:,k:end,:);
        case 5; y = x(:,:,:,:,k-1); z = x(:,:,:,:,k:end);
        otherwise; error('only supported up to 5 dimensions - fix me');
    end
    d = abs(bsxfun(@minus,y,z));
    
    % NaNs don't count to the sum
    if omitnan; d(isnan(d)) = 0; end

    % accumulate sum
    switch dim
        case 1; sumd(k-1,:,:,:,:) = sumd(k-1,:,:,:,:) + sum(d,dim);
                sumd(k:end,:,:,:,:) = sumd(k:end,:,:,:,:) + d;
        case 2; sumd(:,k-1,:,:,:) = sumd(:,k-1,:,:,:) + sum(d,dim);
                sumd(:,k:end,:,:,:) = sumd(:,k:end,:,:,:) + d;
        case 3; sumd(:,:,k-1,:,:) = sumd(:,:,k-1,:,:) + sum(d,dim);
                sumd(:,:,k:end,:,:) = sumd(:,:,k:end,:,:) + d;
        case 4; sumd(:,:,:,k-1,:) = sumd(:,:,:,k-1,:) + sum(d,dim);
                sumd(:,:,:,k:end,:) = sumd(:,:,:,k:end,:) + d;
        case 5; sumd(:,:,:,:,k-1) = sumd(:,:,:,:,k-1) + sum(d,dim);
                sumd(:,:,:,:,k:end) = sumd(:,:,:,:,k:end) + d;
    end

end

% reinsert NaNs so min can ignore them
if omitnan; sumd(isnan(x)) = NaN; end

% indices for min of sum(|d|)
[~,index] = min(sumd,[],dim);

% convert to subs and replace with index
[s1 s2 s3 s4 s5] = ind2sub(szy,1:prod(szy));

switch dim
    case 1; s1 = reshape(index,1,[]);
    case 2; s2 = reshape(index,1,[]);
    case 3; s3 = reshape(index,1,[]);
    case 4; s4 = reshape(index,1,[]);
    case 5; s5 = reshape(index,1,[]);
end

% convert subs back to indicies
index = sub2ind(szx,s1,s2,s3,s4,s5);

% return min elements
y = reshape(x(index),szy);
