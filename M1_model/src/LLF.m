function fLLF = LLF(un, dt, dx, f, lambda, psi2, Nz)
    %% Lax-Friedrichs Fluss nach (2.12)
    % Linker und Rechter Zustand
    ip0c = [1:Nz]'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1];
    ulc = un(:,ip0c); urc = un(:,ip1c);

    % Zentrale Differenz
    fCD = 0.5*(f([urc; psi2(ip1c)])+f([ulc; psi2(ip0c)]));

    % Diffusiver Term
    fLLF = fCD - 0.5 * lambda .* (urc-ulc);
end