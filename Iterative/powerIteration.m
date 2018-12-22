function [lambda, v] = powerIteration(A, v0, maxit)
% This function uses the power method to find the 
% largest eigenvalue and corresponding eigenvector
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to compute eigenvalue & eigenvector of
% v0: initial guess of eigenvector
% maxit: maximum iteration
% OUTPUT:
% lambda: eigenvalue
% v: corresponding eigenvector
%%
v = v0;

for k=1:maxit
    w = A*v;
    v = w / norm(w, 2);
    lambda = v' * A * v;
   
end
end