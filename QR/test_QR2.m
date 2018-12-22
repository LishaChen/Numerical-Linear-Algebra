clear all
close all
clc
%% This is the main function to test QR algorithms for a different matrix
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%% test

% matrix Z
Z = [1,2,3;
    4,5,6;
    7,8,7;
    4,2,3;
    4,2,2];

[m,n] = size(Z);

% Modified GS
[Qm, Rm] = mgs(Z);

% Householder
[Wh,Rh] = house(Z);
Qh = formQ(Wh);
Qh = Qh(:,1:n);
Rh = Rh(1:n,:);

% Matlab QR
[Q,R] = qr(Z);
Q = Q(:,1:n);
R = R(1:n,:);

% Print results
fprintf('Z=\n');
disp(Z);

fprintf('mgs:\n');
fprintf('Q=\n');
disp(Qm);
fprintf('R=\n');
disp(Rm);
fprintf('||A-QR||=%8.2e, ||Q*Q-I|| = %8.2e\n',...
norm(Z-Qm*Rm,2), norm(Qm'*Qm - eye(n)));

fprintf('Householder:\n');
fprintf('Q=\n');
disp(Qh);
fprintf('R=\n');
disp(Rh);
fprintf('||A-QR||=%8.2e, ||Q*Q-I|| = %8.2e\n',...
norm(Z-Qh*Rh,2), norm(Qh'*Qh - eye(n)));

fprintf('Matlab QR: \n');
fprintf('Q=\n');
disp(Q);
fprintf('R=\n');
disp(R)
fprintf('|| A - Q*R || = %8.2e, ||Q*Q-I|| = %8.2e\n',...
    norm(Z-Q*R,2), norm(Q'*Q - eye(n)));

