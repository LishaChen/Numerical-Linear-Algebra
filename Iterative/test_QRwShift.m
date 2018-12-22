clear all
close all
clc
% This function computes eigenvalues by the QR algorithm
% using the shifted algorithm (Wilkinson shift)
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

% QR with shift
% Wilkinson shift
tol = 1E-5;
eigenvalues = QRwShift(A, tol);

err = max(abs(sort(eigenvalues) - sort(eig(A))));
fprintf('eigenvalues:\n')
fprintf('%5.4f\n ', eigenvalues);
fprintf('err %8.5e\n', err);


