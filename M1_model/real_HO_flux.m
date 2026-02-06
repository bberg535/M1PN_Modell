function f = real_HO_flux(uPN, u, dt, dx, method, f, Nz, NPN, source_strength, sigma_a, sigma_s)
% Berechne HO- und LO-Fluss und limitiere durch ein Verfahren aus dem
% Lehrbuch

    % Verschiedene Wellengeschwindigkeiten berechnen
    [lambdaPN,~] = gausslegendre(NPN+1);
    lambdaPN = abs(lambdaPN);
    lambda = 1;

    % Fluss-Berechnung
    LO_psi2 = calc_psi2(u);
    HO_Flux = method(uPN(1:2,:), dt, dx, f, lambdaPN(1:2,:), uPN(3,:), Nz);
    LO_Flux = method(u, dt, dx, f, lambda, LO_psi2, Nz);

    % Limiting
    f = flux_limit(HO_Flux, LO_Flux, Nz, u, LO_psi2, dt, dx , source_strength, sigma_a, sigma_s);
end

