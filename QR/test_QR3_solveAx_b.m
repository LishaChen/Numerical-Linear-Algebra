clear all
close all
clc
%% This is the main function to use QR algorithms to solve Ax = b
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%% test

% form matrix A, vector b
m = 50;
n = 12;
t = (0:m-1)'./(m-1);
A = fliplr(vander(t));
A = A(:, 1:n);

b = cos(4*t);

%% Compute results

% normal equation
xa = (A'*A)\A'*b;

% clgs
[Qc,Rc] = clgs(A);
xb = Rc\Qc'*b;

% mgs
[Qm,Rm] = mgs(A);
xc = Rm\Qm'*b;

% Householder routine
[Wh,Rh] = house(A);
Qh = formQ(Wh);
Qh = Qh(:,1:n);
Rh = Rh(1:n,:);
xd = Rh\Qh'*b;

% Matlab QR
[Q,R] = qr(A);
Q = Q(:,1:n);
R = R(1:n,:);
xe = R\Q'*b;

% Matlab backslash
xf = A\b;

% Matlab SVD
[U, Sigma, V] = svd(A,0);
xg = V/Sigma*U'*b;


%% Print results

fprintf(' Normal                          CLGS                          MGS                         HOUSE\n');
for i=1:n
 fprintf(' %22.15e %22.15e %22.15e %22.15e\n',xa(i),xb(i),xc(i),xd(i))
end;

fprintf(' Matlab QR                          Matlab backslash                          SVD\n');
for i=1:n
fprintf(' %22.15e %22.15e %22.15e\n',xe(i),xf(i),xg(i))
end;
