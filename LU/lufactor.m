function [L,U,P] = lufactor(A)
% This function is the LU factorization with partial pivoting
% PA = LU
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: m*m matrix to be factorized
% OUTPUT:
% L: m*m lower triangular matrix with unit diagonal
% U: m*m upper diagonal matix.
% P: m*m permutation matrix
%%

[m,~] = size(A);

U = A; L=eye(m); P=eye(m);
for k=1:m-1
    [~, idx] = max(abs(U(k:m,k))); % find pivot
    idx = idx + k-1;
    % swap
    uswap = U(idx, k:m);
    U(idx, k:m) = U(k, k:m);
    U(k, k:m) = uswap;
    
    lswap = L(idx, 1:k-1);
    L(idx, 1:k-1) = L(k, 1:k-1);
    L(k, 1:k-1) = lswap;
    
    pswap = P(idx, 1:m);
    P(idx, 1:m) = P(k, 1:m);
    P(k, 1:m) = pswap;
    
    for j = k+1:m
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:m) = U(j,k:m) - L(j,k) * U(k,k:m);
    end
    
end

end