function [x, nit] = gmres0(A, b, x0, maxit, tol)
% This function solves Ax=b using GMRES
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
H = zeros(size(A, 1), 1);
Q = zeros(size(A, 1), 1);
Q(:,1) = b/ norm(b, 2);
nit = maxit;
residOld = norm(b - A*x0, 2);
for n =1:maxit
    v = A*Q(:,n);
    for j = 1:n
        H(j,n) = Q(:,j)' * v;
        v = v - H(j,n) * Q(:,j);        
    end
    H(n+1,n) = norm(v,2);
    Q(:,n+1) = v/ H(n+1, n);
    
    y = (H'*H) \ (norm(b,2)*H'*eye(size(b,1), 1));
    x = Q(:, 1:n) * y;
    resid = norm(H * y - norm(b,2)*eye(size(b,1),1), 2);
    if resid < tol 
        nit = n;
        break
    end
    fprintf('GMRES: n=%3d, || r_n ||_2 = %8.2e, ratio=%8.2e \\\\ \n',...
        n,resid,resid/residOld);
    residOld = resid;
end

end