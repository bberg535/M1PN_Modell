function psi2 = calc_psi2(u)
    % Durch 0 Teilen vermeiden
    eps_rho = 1e-14;
    eps_f = 1e-12;

    % f nach Gleichung (8) berechnen
    f = abs(u(2,:))./(u(1,:) + eps_rho);
    f = min(max(f, 0), 1 - eps_f);

    % Eddington-Tensor in 1D
    eddington = @(fval) ((3 + 4 * fval.^2) ...
        ./ (5 + 2 * sqrt(max(4 - 3 * fval.^2, eps_f))));

    % Zweites Moment berechnen nach Gleichung (5)
    psi2 = eddington(f) .* u(1,:);
end
