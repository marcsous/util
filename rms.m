function out = rms(in,dim)

out = dot(in,in,dim);
out = realsqrt(real(out));