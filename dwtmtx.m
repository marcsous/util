function D = dwtmtx(n,wname)
%DWTMTX Discrete wavelet transform matrix.
%   DWTMTX(N) is the N-by-N matrix of values whose inner
%   product with a column vector of length N yields the
%   discrete wavelet transform of the vector.
%   If X is a column vector of length N, then DWTMTX(N)*X
%   yields the same result as DWT(X); however, DWT(X) is
%   more efficient.
%
%   The inverse discrete wavelet transform matrix is
%   DWTMTX(N)'.

% 'n' input must be a real nonnegative integer scalar
validateattributes(n,{'numeric'},{'real','nonnegative','integer','scalar'},mfilename,'n',1);
if mod(n,2)
    error('n must be even.');
end

% default
if nargin<2
    wname = 'db2';
end

% coefficients
[LO_D HI_D] = wfilters(wname);

% don't entertain oddball situations
if numel(LO_D)>n
    error('filter length (%i) > n is not supported',numel(LO_D));
end

% the matrix
D = zeros(n);

for k = 1:n/2
    D(k,1:numel(LO_D)) = flip(LO_D);
    D(k,:) = circshift(D(k,:),[0 2*k-numel(LO_D)/2-1]);
end
for k = 1+n/2:n
    D(k,1:numel(HI_D)) = flip(HI_D);
    D(k,:) = circshift(D(k,:),[0 2*k-numel(LO_D)/2-1-n]);
end
