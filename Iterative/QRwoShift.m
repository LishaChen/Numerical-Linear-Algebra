function Ak = QRwoShift(A, tol)
% This function computes eigenvalues by the QR algorithm
% using the unshifted algorithm
% Author: Lisha Chen
% Contact: lishachen00@gmail.com
%
% INPUT:
% A: matrix to compute eigenvalue of
% tol: tolerance to stop iteration when eigenvalues converge
% OUTPUT:
% Ak: diagonal are eigenvalues of A
%%
Ak = A;
deltaOld = max(reshape(abs(Ak-diag(diag(Ak))), [], 1));
k=0;
while deltaOld > tol
    k = k+1;
    [Qk, Rk] = qr(Ak);
    Ak = Rk*Qk;
    delta = max(reshape(abs(Ak-diag(diag(Ak))), [],1));
    if mod(k, 10)==0
        fprintf('QR: k=%d : delta=%8.2e, ratio=%5.4f\n',k,delta,delta/deltaOld);
    end
    deltaOld = delta;
end

end