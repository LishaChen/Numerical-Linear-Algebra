clear all
close all
clc
%% This is the main function to test QR algorithms
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%% test

% Form the matrix A to be factorized
m=5;
x = (0:m-1)'./(m-1);
A = fliplr(vander(x));


[Qc, Rc]=clgs(A); % classical GS
[Qm, Rm] = mgs(A); % modified GS
[Wh,Rh] = house(A); % Householder QR
Qh = formQ(Wh); % Householder QR
% Modify Householder QR, change sign of diagnal entries
D = diag(sign(diag(Rh)));
Qh = Qh*D;
Rh = D*Rh;

[Q,R] = qr(A); % matlab QR
% Modify matlab QR, change sign of diagnal entries
D = diag(sign(diag(R)));
Q = Q*D;
R = D*R;

%% compute and print error
fprintf('clgs: m=%d: ||A-QR||=%8.2e, ||Qc-Q|| = %8.2e, ||Rc-R||=%8.2e, ||Qc*Qc-I||=%8.2e \n',...
    m,norm(A-Qc*Rc,2),norm(Qc-Q,2), norm(Rc-R,2),norm(Qc'*Qc-eye(m),2));
fprintf('mgs : m=%d: ||A-QR||=%8.2e, ||Qm-Q|| = %8.2e, ||Rm-R||=%8.2e, ||Qm*Qm-I||=%8.2e \n',...
    m,norm(A-Qm*Rm,2),norm(Qm-Q,2), norm(Rm-R,2),norm(Qm'*Qm-eye(m),2));
fprintf('Householder : m=%d: ||A-QR||=%8.2e, ||Qh-Q|| = %8.2e, ||Rh-R||=%8.2e, ||Qh*Qh-I||=%8.2e \n',...
    m,norm(A-Qh*Rh,2),norm(Qh-Q,2), norm(Rh-R,2),norm(Qh'*Qh-eye(m),2));
fprintf('Matlab QR: || A - Q*R || = %8.2e, ||Q*Q-I|| = %8.2e\n',norm(A-Q*R,2), norm(Q'*Q-eye(m)));


