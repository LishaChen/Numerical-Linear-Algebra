function eigenvalues = QRwShift(A, tol)
% This function computes eigenvalues by the QR algorithm
% using the shifted algorithm (Wilkinson shift)
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to compute eigenvalue of
% tol: tolerance to stop iteration when eigenvalues converge
% OUTPUT:
% eigenvalues: eigenvalues of A
%%
m = size(A, 1);
Ak = A;
deltaOld = max(reshape(abs(Ak-diag(diag(Ak))), [], 1));
k=0;
eigenvalues = [];
dim = m;
while dim >  1
    k = k+1;   
    dt = (Ak(end-1,end) - Ak(end,end))/2;
    mu = Ak(end,end) - (sign(dt)* Ak(end,end-1)^2)/...
        ( abs(dt) + sqrt(dt^2+Ak(end,end-1)^2)); % Wilkinson shift
    
    [Qk, Rk] = qr(Ak - mu * eye(size(Ak, 1)));
    Ak = Rk*Qk + mu * eye(size(Ak, 1));
        
    delta = max(abs([Ak(end, end-1), Ak(end-1, end)]));
    if delta < tol
        eigenvalues = [eigenvalues; Ak(end, end)];
        Ak(:, end) = [];
        Ak(end, :) = [];
        dim = dim-1;
    end

    fprintf('shifted QR: k=%d : shift=%8.2e, dimension=%d, delta=%8.2e, ratio=%5.4f \\\\ \n',...
        k, mu, dim, delta, delta/deltaOld);

    deltaOld = delta;
end
eigenvalues =  [eigenvalues; Ak];
end