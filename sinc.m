function f = sinc(x)

nz = (x~=0);
f(~nz) = 1;
f(nz) = sin(x(nz)) ./ x(nz);
