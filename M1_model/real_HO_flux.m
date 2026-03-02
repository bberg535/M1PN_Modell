function Gstar = real_HO_flux(uPN, u, dt, dx, ilim, flux, nec, NPN, source_strength, sigma_a, sigma_s)
% Berechne HO- und LO-Interface-Fluesse und limitiere den Korrekturfluss.

    nne = nec + 1;
    if size(u, 2) ~= nne || size(uPN, 2) ~= nne
        error('real_HO_flux:dimensionMismatch', ...
            'Expected u and uPN to have %d columns (nec+1).', nne);
    end

    % M1-Fluss mit einheitlicher Wellengeschwindigkeit begrenzen.
    % PN liefert nur die HO-Momente, nicht die Charakteristikgeschwindigkeit.
    lambda = 1;

    % Fluss-Berechnung
    LO_psi2 = calc_psi2(u);
    HO_psi2 = uPN(3, :);

    switch ilim
        case 1
            % MCL: LO=LLF, HO=LW, danach limitierter Korrekturfluss.
            LO_Flux = LLF(u, dt, dx, flux, lambda, LO_psi2, nec);
            HO_Flux = LW(uPN(1:2, :), dt, dx, flux, lambda, HO_psi2, nec);
            Gstar = flux_limit(HO_Flux, LO_Flux, nec, u, LO_psi2, dt, dx, source_strength, sigma_a, sigma_s);
        case -1
            % Reiner LO-Flux.
            Gstar = LLF(u, dt, dx, flux, lambda, LO_psi2, nec);
        otherwise
            % Reiner HO-Target-Flux auf Basis der PN-Momente.
            Gstar = LW(uPN(1:2, :), dt, dx, flux, lambda, HO_psi2, nec);
    end
end

