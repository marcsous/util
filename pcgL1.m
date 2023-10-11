function [x lambda resvec] = pcgL1(A,b,sparsity,tol,maxit,Q)
% [x lambda resvec] = pcgL1(A,b,sparsity,tol,maxit,Q)
%
% Solves the following problem via ADMM:
%
%   minimize (1/2)||Ax-b||_2^2 + lambda*||Qx||_1
%
% -A is an n x n symmetric matrix (function handle)
% -sparsity is the fraction of zeros (0.2=20% zeros)
% -tol/maxit are tolerance and max. no. iterations
% -Q is a wavelet transform Q*x and Q'*x (e.g. HWT.m)
%
% -lambda that yields the required sparsity (scalar)
% -resvec is the residual at each iteration (vector)
%
% Derived from lasso_lsqr.m:
% https://web.stanford.edu/~boyd/papers/admm/lasso/lasso_lsqr.html
%
%% check arguments
if nargin<3 || nargin>6
    error('Wrong number of input arguments');
end
if nargin<4 || isempty(tol)
    tol = 1e-6;
end
if nargin<5 || isempty(maxit)
    maxit = 100;
end
if nargin<6 || isempty(Q)
    Q = 1;
end
if isnumeric(A)
    A = @(arg) A * arg;
end
if ~iscolumn(b)
    error('b must be a column vector');
end
if ~isscalar(sparsity) || sparsity<0 || sparsity>1
    error('sparsity must be a scalar between 0 and 1.');
end

% check A is [n x n]
n = numel(b);
try
    Ab = A(b);
catch
    error('A*x failed when passed a vector of length %i.',n);
end   
if ~iscolumn(Ab) || numel(Ab)~=n
    error('A*x did not return a vector of length %i.',n);
end

% check positive definiteness (50% chance of catching it)
bAb = b'*Ab;
if abs(imag(bAb)) > sqrt(n)*eps(bAb) || real(bAb) < -eps(bAb)
    warning('Matrix operator A may not be positive definite.');
end

%% solve
t = tic();

alpha = (b'*b) / bAb;
x = alpha * b;

z = zeros(size(b),'like',b);
u = zeros(size(b),'like',b);

for iter = 1:maxit 
 
    % x-update
    Ak = @(x)A(x) + x * (iter>1);
    bk = b + (z - u);
    [x flag] = minres(Ak,bk,[],[],[],[],x);
    
    if flag>4
        warning('minres returned error flag %i.',flag);
    end
    
    % z-update
    z = Q * (x + u);
    
    [z lambda(iter)] = shrinkage(z, sparsity);
    
    z = Q' * z;
    
    u = u + (x - z);
    
    % check convergence
    resvec(iter) = norm(x-z) / norm(x);

    if resvec(iter) < tol
        break;
    end
    
end
x = z;

% report convergence
if nargout<2
    if iter < maxit
        fprintf('%s converged at iteration %i to a solution with relative residual %.1e. ',mfilename,iter,resvec(iter));
    else
        fprintf('%s stopped at iteration %i with relative residual %.1e without converging to the desired tolerance %.1e. ',mfilename,iter,resvec(iter),tol);
    end
    toc(t);
end

%% shrink based on sparsity, return lambda
function [z lambda] = shrinkage(z,sparsity)

absz = abs(z);

v = sort(absz,'ascend');

lambda = interp1(v,sparsity*numel(z),'linear',0);

z = sign(z) .* max(absz-lambda, 0); % complex ok
