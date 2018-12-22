function Q = formQ(W)
% This function forms Q from Householder vectors
% from the function hessenberg
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% W: mxn lower triangular matrix, 
%      columns are Householder vectors
% OUTPUT:
% Q: mxm unitary matrix
%%

[m,n]=size(W);
Q = eye(m);

for k=n:-1:1
    Q(k:m,:) = Q(k:m,:) - 2*W(k:m,k)*((W(k:m,k))'*Q(k:m,:));
end

end