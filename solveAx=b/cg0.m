function [x,nit] = cg0( A,b,x0, maxit, tol )
% This function solves Ax=b using Conjugate Gradient algorithm
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: nonsingular matrix
% b: right-hand-side
% x0: initial guess
% maxit: maximum iteration
% tol: residual convergence tolerance
% OUTPUT:
% x: solution
% nit: number of iterations
%%
x = zeros(size(b));
r = b;
p = r;
residOld = norm(b - A*x0, 2);

for n=1:maxit
    alpha = (r' * r) / (p'*A*p);
    x = x + alpha*p;
    rn= r - alpha*A*p;
    beta = (rn' * rn) / (r'*r);
    r = rn;
    p = r + beta*p;
    resid = norm(b - A*x, 2);
    if resid < tol 
        nit = n;
        break
    end
    fprintf('CG: n=%3d, || r\\_n ||\\_2 = %8.2e, ratio=%8.2e \\\\ \n',...
        n,resid,resid/residOld);
    residOld = resid;
    
end

end