clear all
close all
clc

% This function computes an implicit  representation
% of the similarity transform to an upper Henssenberg matrix H
% A = QHQ*
% Author: Lisha Chen
% Contact: lishachen00@gmail.com

%% test
m=5;
A = zeros(m);

for i=1:m
    for j=1:m
        
        if i==j
            A(i,j) = 9;
        else
            A(i,j) = 1/(i+j);
        end
     end
end

 [W,H] = hessenberg(A)
 Q = formQ(W)
 
 % compute error
 norm(Q'*Q - eye(m), 2)
 norm(A - Q*H*Q', 2)
 
 


