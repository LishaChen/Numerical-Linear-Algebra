function [W,H] = hessenberg(A)
% This function computes an implicit  representation
% of the similarity transform to an upper Henssenberg matrix H
% A = QHQ*
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: m*m matrix
% OUTPUT:
% W: m*m lower triangular matrix with columns the Housholder vectors
% H: m*m upper Hessenberg matrix
%%

[m,n] = size(A);
W = zeros(m,n);
for k=1:m-2
    x = A(k+1:m,k);
    W(k+1:m,k) = sign(x(1)) * norm(x, 2) * eye(size(x,1),1) + x;
    W(k+1:m,k) = W(k+1:m,k) ./ norm(W(k+1:m,k), 2);
    A(k+1:m,k:m) = A(k+1:m,k:m) - 2*W(k+1:m,k)*(W(k+1:m,k)' * A(k+1:m,k:m));
    A(1:m,k+1:m) = A(1:m,k+1:m) - 2* A(1:m,k+1:m)* W(k+1:m,k)*W(k+1:m,k)' ;
end

H = A;

end