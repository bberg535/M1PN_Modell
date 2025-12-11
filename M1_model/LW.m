function u = LW(un, dt, dx, f, lambda, psi2, Nz)
    %% Lax-Wendroff Fluss nach (2.20)
    % Linker und Rechter Zustand
    ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1];
    ulc = un(:,ip0c); urc = un(:,ip1c);

    % Zentrale Differenz
    fCD = 0.5*(f([urc; psi2(ip1c)])+f([ulc; psi2(ip0c)]));

    % Diffusiver Term
    u = fCD - 0.5 * dt/dx * (lambda.^2) .* (urc-ulc);
end