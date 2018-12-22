function [W, R] = house(A)
% This function implements the Householder reflector:
%         Implicitly Compute the reduced QR factorization
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: mxn matrix
% OUTPUT:
% W: mxn lower triangular matrix, 
%       columns are Householder vectors
% R: nxn matrix upper triangular
%%
[m,n] = size(A);

W = zeros(m,n);

for k=1:n
    x = A(k:m,k);
    W(k:m,k) = sign(x(1))*norm(x,2)*eye(size(x,1),1)+x;
    W(k:m,k) = W(k:m,k) ./ norm(W(k:m,k), 2);
    A(k:m,k:n) = A(k:m,k:n) - 2*W(k:m,k)*((W(k:m,k))'*A(k:m,k:n));
    
end

R = A;

end