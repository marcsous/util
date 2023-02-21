function [w] = gausswin(M,alpha)
if nargin<2; alpha = 2.5; end
n = -(M-1)/2 : (M-1)/2;
w = exp((-1/2) * (alpha * n/((M-1)/2)) .^ 2)';