function [lambda, v] = rayleighQuotientIteration(A, lambda0, v0, maxit)
% This function uses the Rayleigh Quotient method to find an
% eigenvalue and corresponding eigenvector
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to compute eigenvalue & eigenvector of
% lambda0: initial guess of eigenvalue
% v0: initial guess of eigenvector
% maxit: maximum iteration
% OUTPUT:
% lambda: eigenvalue
% v: corresponding eigenvector
%%

m = size(A, 1);

% Rayleigh Quotient method
v = v0 / norm(v0, 2);
lambda = lambda0;
% lambda = v' * A * v; 
% can also use RQ to compute initial guess of eigenvalue

for k=1:maxit
    w =   (A- lambda*eye(m))\ v; 
    v = w / norm(w,2);
    lambda = v' * A * v;
    
end
end