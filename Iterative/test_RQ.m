clear all
close all
clc

% This function uses the Rayleigh Quotient method to find an
% eigenvalue and corresponding eigenvector
% Author: Lisha Chen
% Contact: lishachen00@gmail.com

%% test
% dimension m
m = 10;

% matrix A
A = zeros(m,m);
for i=1:m
    if i > 1
        A(i,i-1) = -1;
    end
    A(i,i) = 4+i;
    if i<m
        A(i,i+1) = -1;
    end
end


% matlab function to compute eigenvalue and eigenvector
[V,D] = eig(A);
lambda1 = D(6,6); % matlab computed lambda1
v1 = V(:,end); % matlab computed q1

% Rayleigh Quotient method
% initial v and lambda
v0 = ones(m, 1);
lambda0 = 10.5;
% max iteration
maxit = 5;
[lambda, v] = rayleighQuotientIteration(A, lambda0, v0, maxit);

