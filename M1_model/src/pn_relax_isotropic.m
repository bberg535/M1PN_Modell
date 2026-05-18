function u_out = pn_relax_isotropic(u_in, N, dt, sigma_a, sigma_s, q)
% Kollision/Absorption und Quelle fuer das PN-Modell
    Nz = size(u_in, 2);
    sigma_a = reshape(sigma_a, 1, []);
    sigma_s = reshape(sigma_s, 1, []);
    if isscalar(q)
        q = q * ones(1, Nz);
    else
        q = reshape(q, 1, []);
    end
    if numel(sigma_a) ~= Nz || numel(sigma_s) ~= Nz || numel(q) ~= Nz
        error('pn_relax_isotropic:dimensionMismatch', ...
            'sigma_a, sigma_s and q must have length %d.', Nz);
    end

    % Transformations- und Basisvektoren/Matrizen
    [~, w, R, L] = pn_basis(N);

    % Diagonalisieren
    v = L * u_in;

    % Berechne rho, um Reaktions- und Quellterme hinzuzufuegen
    rho = (w.' * v);

    % rho^{n+1} aus der Momentgleichung (implizit in sigma_a)
    rho_new = (rho + dt * q) ./ (1 + dt * sigma_a);

    % denominators
    denom = 1 + dt * (sigma_a + sigma_s);

    % isotrope Beitraege in jeder Richtung: (sigma_s / 2) * rho_new + q / 2
    iso = 0.5 * (sigma_s .* rho_new + q);

    % update in Winkelraum (implizit in sigma_t)
    v_new = (v + dt * (ones(N + 1, 1) * iso)) ./ (ones(N + 1, 1) * denom);

    % zurueck in Momentraum
    u_out = R * v_new;

    if all(sigma_a == 0) && all(q == 0)
        tol = 1.0e-11 * max(1.0, max(abs(u_in(:))));
        rho_err = max(abs(u_out(1, :) - u_in(1, :)));
        if rho_err > tol
            error('pn_relax_isotropic:rhoNotConserved', ...
                ['Isotropic scattering must leave u_0 unchanged for sigma_a = 0, q = 0. ' ...
                 'Observed max |delta u_0| = %.6e.'], rho_err);
        end

        expected = u_in;
        expected(2:end, :) = u_in(2:end, :) ./ (ones(N, 1) * denom);
        higher_err = max(abs(u_out(2:end, :) - expected(2:end, :)), [], 'all');
        if higher_err > tol
            error('pn_relax_isotropic:higherMomentsMismatch', ...
                ['Isotropic scattering must damp only l >= 1 moments by 1/(1 + dt*sigma_s). ' ...
                 'Observed max mismatch = %.6e.'], higher_err);
        end
    end
end
