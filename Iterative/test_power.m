clear all
close all
clc

% This function uses the power method to find the 
% largest eigenvalue and corresponding eigenvector
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
lambda1 = D(end); % matlab computed lambda1
v1 = V(:,end); % matlab computed q1

% power method
v0 = ones(m, 1); % initial v
maxit = 25; % max iteration
[lambda, v] = powerIteration(A, v0, maxit);




