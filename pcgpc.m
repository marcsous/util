function [x,flag,relres,iter,resvec] = pcgpc(A,b,tol,maxit,M1,M2,x0)
% Preconditioned conjugate gradient (pcg) with modifications to
% support the use of a penalty on Im(x) aka phase constraint.
%
% Meant to be used with anonymous functions, A = @(x)myfunc(x),
% where myfunc(x) returns:
% - A*x           : to minimize ||b-Ax||^2
% - A*x+λ*x       : to minimize ||b-Ax||^2 + λ^2||x||^2
% - A*x+λ*i*Im(x) : to minimize ||b-Ax||^2 + λ^2||Im(x)||^2
%
% References:
% - An Introduction to the CG Method Without the Agonizing Pain
%   (Jonathan Richard Shewchuk, Carnegie Mellon University 1994)
% - Partial Fourier Partially Parallel Imaging
%   (Mark Bydder and Matthew D. Robson, MRM 2005; 53: 1393)
%
% Modifications:
% - uses the real part only of the dot products
% - allows multiple RHS vectors (b = [b1 b2 ...])
% - mostly compatible with pcg (except M2 and flag)
%
% Usage: see Matlab's pcg function (help pcg)

% check arguments
if nargin<2; error('Not enough input arguments.'); end
if ~ismatrix(b); error('b argument must be a column vector or 2d array'); end
if ~exist('tol') || isempty(tol); tol = 1e-6; else; tol = gather(tol); end
if ~exist('maxit') || isempty(maxit); maxit = min(20,size(b,1)); end
if ~exist('M1') || isempty(M1); M1 = @(arg) arg; end
if exist('M2') && ~isempty(M2); error('M2 argument not supported'); end
validateattributes(tol,{'numeric'},{'scalar','nonnegative','finite'},'','tol');
validateattributes(maxit,{'numeric'},{'scalar','nonnegative','integer'},'','maxit');

% not intended for matrix inputs but they are supported
if isnumeric(A); A = @(arg) A * arg; end
if isnumeric(M1); M1 = @(arg) M1 \ arg; end

% initialize
iter = 1;
flag = 1;
if ~exist('x0') || isempty(x0)
    r = b;
    x = zeros(size(b),'like',b);
else
    if ~isequal(size(x0),size(b))
        error('x0 must have length %i to match the problem size.',size(b,1));
    end
    x = x0;
    r = A(x);   
    if ~isequal(size(r),size(b))
        error('A(x) must return a column vector of length %i to match the problem size.',numel(b));
    end  
    r = b - r;
end
d = M1(r);
if ~isequal(size(d),size(b))
    error('M1(x) must return a column vector of length %i to match the problem size.',numel(b));
end
delta0 = vecnorm(b);
delta_new = real(dot(r,d));
resvec(iter,:) = vecnorm(r);

% min norm solution
xmin = x;
imin = zeros(1,size(b,2));

% main loop
while maxit
    
    iter = iter+1;
	clear q; q = A(d);
    if ~isequal(size(b),size(q))
        error('b is [%s] but matrix vector product is [%s].',num2str(size(b)),num2str(size(q)));
    end
	alpha = delta_new./real(dot(d,q));
    
    % unsuccessful termination
    if ~all(isfinite(alpha)); flag = 4; break; end
    
	x = x + alpha.*d;
    r = r - alpha.*q;
    
    % recalculate residual occasionally?
    if 0 %mod(iter,50)==0
        r = b - A(x);
    end

    % residual vectors
    resvec(iter,:) = vecnorm(r);

    % keep best solution
    ok = resvec(iter,:) < min(resvec(1:iter-1,:));
    if any(ok)
        xmin(:,ok) = x(:,ok);
        imin(:,ok) = iter;
    end
    
    % successful termination
    if all(resvec(iter,:)<tol*delta0); flag = 0; break; end

    % unsuccessful termination
    if iter>maxit; flag = 1; break; end

    clear q; q = M1(r);
	delta_old = delta_new;
    delta_new = real(dot(r,q));
   
    % unsuccessful termination
    if all(delta_new<=0); flag = 4; break; end

    beta = delta_new./delta_old;
    d = q + beta.*d;

end

% min norm solution
x = xmin;
iter = imin;

% min residual norm
relres = resvec(iter,:)./delta0;

% only display if flag not supplied
if nargout<2
    for k = 1:size(b,2)
        fprintf('pcg terminated at iteration %i (flag %i): relres = %e.\n',imin(k),flag,relres(k)); 
    end
end
