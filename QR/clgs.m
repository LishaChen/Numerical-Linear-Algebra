function [Qc, Rc]=clgs(A)
% This function is the classical Gram-Schmidt algorithm
% to compute the reduced QR factorization: A=QR
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to be factorized
% OUTPUT:
% Qc: m*n matrix with orthonormal columns
% Rc: n*n matrix, upper triangular
%%
[m,n] = size(A);
Rc = zeros(n,n);
Qc = zeros(m,n);
I=1:m;
for j = 1:n
    vj = A(I,j);
    for i=1:j-1
        Rc(i,j) = dot(Qc(I,i), A(I,j));
        vj=vj-Rc(i,j)*Qc(I,i);
    end
Rc(j,j) = norm(vj, 2);
Qc(I,j) = vj/Rc(j,j);
end

end