function psi2 = calc_psi2(u)
    % Durch 0 Teilen vermeiden
    eps_rho = 1e-14;

    % f nach Gleichung (8) berechnen
    f   = abs(u(2,:))./(u(1,:) + eps_rho);

    % Eddington-Tensor in 1D
    eddington = @(f) ((3 + 4 * f.^2)./(5 + 2*sqrt(4-3*f.^2)));

    % Zweites Moment berechnen nach Gleichung (5)
    psi2 = eddington(f).*u(1,:);
end