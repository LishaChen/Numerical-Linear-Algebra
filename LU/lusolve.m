function x = lusolve(b,L,U,P)
% This function is to solve Ax = b using the results from LU factor
% PA = LU
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT: 
% L,U,P that PA=LU
% b: m*1 vector, Ax=b
% OUTPUT: 
% x: m*1 solution to Ax=b
%%
[m,n] = size(L);

Ux = P*b;
for i =2:m
    res = L(i,1:i-1)*Ux(1:i-1);
    Ux(i) = Ux(i) - res;
end

x = Ux;
x(m) = x(m)/U(m,m);

for i =m-1:-1:1
    res = U(i,i+1:m)*x(i+1:m);
    x(i) = (x(i) - res)/U(i,i);
end


end