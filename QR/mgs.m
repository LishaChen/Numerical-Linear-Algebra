function [Qm, Rm] = mgs(A)
% This function is the modified Gram-Schmidt algorithm
% to compute the reduced QR factorization: A=QR
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to be factorized
% OUTPUT:
% Qm: m*n matrix with orthonormal columns
% Rm: n*n matrix, upper triangular
%%
[m,n] = size(A);
Rm = zeros(n,n);
Qm = zeros(m,n);

v = A;

for i = 1:n
    Rm(i,i) = norm(v(:,i), 2);
    Qm(:,i) = v(:,i)/Rm(i,i);
    
    for j=i+1:n
        Rm(i,j) = dot(Qm(:,i), v(:,j));
        v(:,j)=v(:,j)-Rm(i,j)*Qm(:,i);
    end

end

end