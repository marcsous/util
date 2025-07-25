function out = rms(in,dim)

if nargin==1
    dim = find(size(in)>1,1,'first');
end

out = real(dot(in,in,dim));
out = out/size(in,dim);
out = realsqrt(out);
