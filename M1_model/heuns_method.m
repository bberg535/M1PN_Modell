function [u, uPN] = heuns_method(u, uPN, dt, dz, flux, nec, sigma_a, sigma_s, ilim, N, source_strength, ~)
    nne = nec + 1;
    if size(u, 2) ~= nne || size(uPN, 2) ~= nne
        error('heuns_method:dimensionMismatch', ...
            'Expected u and uPN to have %d columns (nec+1).', nne);
    end

    u = enforce_periodic(u, nec);
    uPN = enforce_periodic(uPN, nec);

    %% Stufe 1
    rhs1 = build_rhs(uPN, u, dt, dz, ilim, flux, nec, N, source_strength, sigma_a, sigma_s);
    u1 = u + dt * rhs1;
    u1 = enforce_periodic(u1, nec);
    u1 = enforce_realizable(u1, nec);

    k1PN = advance_pn(uPN, N, ilim, dt, dz);
    uPN1 = uPN + dt * k1PN;
    uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, source_strength);
    uPN1 = sync_pn_to_m1(uPN1, u1, nec);
    uPN1 = enforce_periodic(uPN1, nec);

    %% Stufe 2
    rhs2 = build_rhs(uPN1, u1, dt, dz, ilim, flux, nec, N, source_strength, sigma_a, sigma_s);
    u2 = u1 + dt * rhs2;
    u2 = enforce_periodic(u2, nec);
    u2 = enforce_realizable(u2, nec);

    k2PN = advance_pn(uPN1, N, ilim, dt, dz);
    uPN2 = uPN1 + dt * k2PN;
    uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, source_strength);
    uPN2 = sync_pn_to_m1(uPN2, u2, nec);
    uPN2 = enforce_periodic(uPN2, nec);

    %% Zeitschritt
    u = 0.5 * (u + u2);
    u = enforce_periodic(u, nec);
    u = enforce_realizable(u, nec);

    uPN = 0.5 * (uPN + uPN2);
    uPN = sync_pn_to_m1(uPN, u, nec);
    uPN = enforce_periodic(uPN, nec);

    if any(~isfinite(u(:))) || any(~isreal(u(:)))
        error('heuns_method:invalidState', ...
            'M1 update produced non-finite or complex values.');
    end
end

function rhs = build_rhs(uPN, u, dt, dz, ilim, flux, nec, NPN, source_strength, sigma_a, sigma_s)
    im1c = [nec, 1:nec-1];
    Gstar = real_HO_flux(uPN, u, dt, dz, ilim, flux, nec, NPN, source_strength, sigma_a, sigma_s);

    rhs = zeros(2, nec + 1);
    source = source_term(nec, source_strength);
    react = reaction(u(:, 1:nec), sigma_a(1:nec), sigma_s(1:nec));
    rhs(:, 1:nec) = -(Gstar(:, 1:nec) - Gstar(:, im1c)) / dz + source - react;
    rhs(:, nec + 1) = rhs(:, 1);
end

function u = enforce_periodic(u, nec)
    u(:, nec + 1) = u(:, 1);
end

function uPN = sync_pn_to_m1(uPN, u, nec)
    uPN(1:2, 1:nec) = u(:, 1:nec);
    uPN(1:2, nec + 1) = u(:, 1);
end

function u = enforce_realizable(u, nec)
    eps_rho = 0;
    eps_f = 1e-12;

    rho = u(1, 1:nec);
    mom = u(2, 1:nec);

    rho = max(rho, eps_rho);
    mom_lim = (1 - eps_f) .* rho;
    mom = min(max(mom, -mom_lim), mom_lim);

    u(1, 1:nec) = rho;
    u(2, 1:nec) = mom;
    u = enforce_periodic(u, nec);
end
