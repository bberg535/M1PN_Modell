function u_out = pn_relax_isotropic(u_in, N, dt, sigma_a, sigma_s, q)
% Kollision/Absorption und Quelle für das PN-Modell
    Nz = size(u_in,2);

    % Transformations- und Basisvektoren/Matrizen
    [xi,w] = gausslegendre(N+1);
    P = legpoly_eval(xi, N);
    G = diag(2./(2.*(0:N)+1));
    R = P*diag(w);
    L = (G \ P).';

    % Diagonalisieren
    v = L * u_in;

    % Berechne rho, um Reaktions- und Quellterme hinzuzufügen
    rho = (w.' * v);

    % rho^{n+1} aus der Momentgleichung (implizit in sigma_a)
    rho_new = (rho + dt*q) ./ (1 + dt*sigma_a);   % 1×Nz

    % denominators
    denom = 1 + dt*(sigma_a + sigma_s);           % 1×Nz

    % isotrope Beiträge in jeder Richtung:
    % (sigma_s/2)*rho_new + (q/2)
    iso = 0.5*(sigma_s .* rho_new + q);           % 1×Nz

    % update in Winkelraum (implizit in sigma_t)
    v_new = (v + dt * (ones(N+1,1) * iso)) ./ (ones(N+1,1) * denom);   % K×Nz

    % zurück in Momentraum
    u_out = R * v_new;                            % (N+1)×Nz
end
