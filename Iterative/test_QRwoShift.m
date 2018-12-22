clear all
close all
clc
% This function computes eigenvalues by the QR algorithm
% using the unshifted algorithm
% Author: Lisha Chen
% Contact: lishachen00@gmail.com

%% test
% dimension m
m = 11;

% matrix A
A = zeros(m,m);
for i=1:m
    if i > 1
        A(i,i-1) = -1;
    end
    A(i,i) = 2;
    if i<m
        A(i,i+1) = -1;
    end
end

% QR without shift

tol = 1E-5;
Ak = QRwoShift(A, tol);

eigenvalues = diag(Ak);
err = max(abs(sort(diag(Ak)) - sort(eig(A))));
fprintf('eigenvalues:\n')
fprintf('%5.4f\n ', eigenvalues);
fprintf('err %8.5e\n', err);


