function Gstar = flux_limit(HO_Flux, LO_Flux, nec, u, psi2, ~, ~, ~, ~, ~)
    ip0c = 1:nec;
    ip1c = [2:nec, 1];
    im1c = [nec, 1:nec-1];

    if size(HO_Flux, 2) ~= nec || size(LO_Flux, 2) ~= nec
        error('flux_limit:dimensionMismatch', ...
            'Expected HO/LO interface fluxes with nec=%d columns.', nec);
    end

    u_act = u(:, 1:nec);
    psi2_act = psi2(1:nec);
    flux = [u_act(2, :); psi2_act];

    % Antidiffusiver Zielfluss (MCL-konform: LO minus HO)
    fAe = LO_Flux - HO_Flux;
    % Die erste M1-Gleichung (rho_t + (u1)_x = ...) bleibt im LO-M1-Flux.
    % PN-Korrektur wird nur auf den Impulsfluss (zweite Komponente) angewandt.
    fAe(1, :) = 0;

    % Lokale Schranken auf Knotenebene
    umax = max(u_act(:, im1c), max(u_act(:, ip0c), u_act(:, ip1c)));
    umin = min(u_act(:, im1c), min(u_act(:, ip0c), u_act(:, ip1c)));

    % Low-order Bar-States auf Interfaces
    wbar = 0.5 * (u_act(:, ip0c) + u_act(:, ip1c)) ...
         - 0.5 * (flux(:, ip1c) - flux(:, ip0c));

    % Limiter fuer den antidiffusiven Fluss
    fAe = min(max(0, fAe), min(umax(:, ip0c) - wbar, wbar - umin(:, ip1c))) ...
        + max(min(0, fAe), max(umin(:, ip0c) - wbar, wbar - umax(:, ip1c)));

    % Vollstaendiger limitierter Interface-Fluss
    Gstar = LO_Flux - fAe;
end
