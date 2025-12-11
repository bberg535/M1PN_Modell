function uPN_new = ls_fit_PN_legendre(uPN_HO, u_M1)
%LS_FIT_PN_LEGENDRE  Least-Squares-Fitting für PN (Legendre-Basis)
%
%   uPN_new = ls_fit_PN_legendre(uPN_HO, u_M1)
%
%   Eingaben:
%       uPN_HO : K×Nz   PN-Koeffizienten nach dem "High-Order"-Schritt
%                        (z.B. nach Heun / SIU / PN-Update)
%                 Annahme: Legendre-Basis mit
%                     uPN_HO(1,:) = u0  ~ U0  (Dichte)
%                     uPN_HO(2,:) = u1  ~ U1  (Fluss)
%
%       u_M1   : 2×Nz   M1-Momente nach dem M1-Schritt:
%                     u_M1(1,:) = U0_new  (Dichte)
%                     u_M1(2,:) = U1_new  (Fluss)
%
%   Ausgabe:
%       uPN_new: K×Nz   gefittete PN-Koeffizienten, so dass
%                 - die ersten beiden Momente zu u_M1 passen:
%                       uPN_new(1,:) = u_M1(1,:)
%                       uPN_new(2,:) = u_M1(2,:)
%                 - uPN_new möglichst nah an uPN_HO liegt
%                   im Sinne von ||uPN_new - uPN_HO||^2 (Least Squares).
%
%   Mathematisch:
%       Für jede Zelle j lösen wir
%
%           min_v ||v - v_HO||^2
%           s.t.  M v = U_new,
%
%       mit  M = [1 0 0 ...;
%                 0 1 0 ...]   (2×K),
%       v_HO  = uPN_HO(:,j),
%       U_new = u_M1(:,j).
%
%   Implementiert wird die geschlossene Lösung:
%
%       v_new = v_HO + M' (M M')^{-1} (U_new - M v_HO)
%
%   ohne lsqlin, also ohne Optimization Toolbox.

    % --- Dimensionen prüfen ------------------------------------------------
    [K, Nz] = size(uPN_HO);
    if size(u_M1, 1) ~= 2
        error('u_M1 muss die Größe 2×Nz haben (erste Zeile U0, zweite Zeile U1).');
    end
    if size(u_M1, 2) ~= Nz
        error('Anzahl Zellen in uPN_HO und u_M1 muss übereinstimmen.');
    end

    % --- Momentmatrix M für Legendre-PN: U0 = u0, U1 = u1 ------------------
    % M ist 2×K:
    %   [U0; U1] = M * v  mit  v = [u0; u1; u2; ...; u_{K-1}]
    M = zeros(2, K);
    M(1,1) = 1;   % U0 = u0
    M(2,2) = 1;   % U1 = u1

    % (M M^T)^{-1} ist 2×2
    MMt_inv = inv(M * M.');

    % --- Ausgabeinitialisierung -------------------------------------------
    uPN_new = zeros(K, Nz);

    % --- Zellweise Least-Squares-Projektion --------------------------------
    for j = 1:Nz
        vHO  = uPN_HO(:, j);   % K×1 : PN-Koeffizienten nach HO-Schritt
        Utar = u_M1(:, j);     % 2×1 : gewünschte Momente [U0; U1] aus M1

        % Korrektur = M' (MM')^{-1} (Utar - M vHO)
        correction   = M' * (MMt_inv * (Utar - M * vHO));
        uPN_new(:,j) = vHO + correction;
    end
end
