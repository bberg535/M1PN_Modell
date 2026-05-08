function fLW = LW(un, dt, dx, f, ~, psi2, Nz)
    %% Lax-Wendroff/Richtmyer-Fluss fuer das 2x2-Momentensystem
    % Linker und rechter Zustand an den Interfaces
    ip0c = (1:Nz)';
    ip1c = [2:Nz 1]';
    ulc = un(:, ip0c);
    urc = un(:, ip1c);

    flux_l = f([ulc; psi2(ip0c)]);
    flux_r = f([urc; psi2(ip1c)]);
    u_half = 0.5 * (ulc + urc) - 0.5 * dt / dx * (flux_r - flux_l);

    % Wenn psi2 aus der M1-Schliessung stammt, wird die Halbzeit-Flussauswertung
    % konsistent aus dem vorhergesagten Zustand rekonstruiert. Fuer extern
    % gelieferte High-Order-Momente (z.B. PN) bleibt psi2 ein Hilfsfeld.
    psi2_closure = calc_psi2(un);
    closure_mismatch = max(abs(psi2(:) - psi2_closure(:)));
    if closure_mismatch <= 1e-10 * max(1, max(abs(psi2_closure(:))))
        %u_half = enforce_realizable_half_step(u_half);
        psi2_half = calc_psi2(u_half);
    else
        psi2_half = 0.5 * (psi2(ip0c) + psi2(ip1c));
    end

    fLW = f([u_half; psi2_half]);
end

function u = enforce_realizable_half_step(u)
    rho = max(u(1, :), 0.0);
    mom_lim = (1 - 1e-12) .* rho;
    mom = min(max(u(2, :), -mom_lim), mom_lim);
    u(1, :) = rho;
    u(2, :) = mom;
end
