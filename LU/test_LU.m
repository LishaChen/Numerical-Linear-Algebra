clear all
close all
clc
% This function is to solve Ax = b using the results from LU factor
% PA = LU
% Author: Lisha Chen
% Contact: lishachen00@gmail.com

%% test
A = [2,1,1,0;
        4,3,3,1;
        8,7,9,5;
        6,7,9,8];
b = [7;23;69;79];

[L,U,P] = lufactor(A) % LU factorization with permutation
norm(P*A-L*U) % error of PA=LU

x = lusolve(b,L,U,P) % solution x to Ax=b
norm(A*x-b) % error of Ax=b

