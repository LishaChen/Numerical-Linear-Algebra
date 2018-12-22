clear all
close all
clc
% This function solves Ax=b using GMRES
% Author: Lisha Chen
% Contact: lishachen00@gmail.com

%% test
% dimension m
m = 100;

% matrix A
A = zeros(m,m);
for i=1:m
    if i > 1
        A(i,i-1) = 1;
    end
    A(i,i) = 4;
    if i<m
        A(i,i+1) = 2;
    end
end
% vector b
b = zeros(m,1);
for i =1:m
    b(i) = 1+i/m;
end

% parameters
x0 = b;
tol = 1e-5;
maxit = 30;

[x, nit] = gmres0(A, b, x0, maxit, tol);

xTrue = A\b;

fprintf('FINAL: Solution from GMRES: (nit=%d)\n',nit);
fprintf('2-norm error=%8.2e\n',norm(x-xTrue,2));
fprintf('max-norm error=%8.2e\n',norm(x-xTrue,inf));
