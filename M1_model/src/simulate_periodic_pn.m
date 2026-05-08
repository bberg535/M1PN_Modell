function uPN = simulate_periodic_pn(rho0, dt, dz, num_steps, N, ilim, sigma_a, sigma_s, q)
%SIMULATE_PERIODIC_PN Advance the periodic PN model with SSP Heun time stepping.

    nec = numel(rho0);
    nne = nec + 1;

    sigma_a = extend_periodic_field(sigma_a, nec);
    sigma_s = extend_periodic_field(sigma_s, nec);
    q = extend_periodic_field(q, nec);

    uPN = zeros(N + 1, nne);
    uPN(1, 1:nec) = reshape(rho0, 1, []);
    uPN(:, nne) = uPN(:, 1);

    for step = 1:num_steps
        k1 = advance_pn(uPN, N, ilim, dt, dz);
        uPN1 = uPN + dt * k1;
        uPN1(:, nne) = uPN1(:, 1);
        uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, q);
        uPN1(:, nne) = uPN1(:, 1);

        k2 = advance_pn(uPN1, N, ilim, dt, dz);
        uPN2 = uPN1 + dt * k2;
        uPN2(:, nne) = uPN2(:, 1);
        uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, q);
        uPN2(:, nne) = uPN2(:, 1);

        uPN = 0.5 * (uPN + uPN2);
        uPN(:, nne) = uPN(:, 1);

        if any(~isfinite(uPN(:))) || any(~isreal(uPN(:)))
            error('simulate_periodic_pn:invalidState', ...
                'PN state became invalid at step %d.', step);
        end
    end
end

function field = extend_periodic_field(field, nec)
    if isscalar(field)
        field = field * ones(1, nec);
    else
        field = reshape(field, 1, []);
    end

    if numel(field) ~= nec
        error('simulate_periodic_pn:dimensionMismatch', ...
            'Expected field length %d.', nec);
    end

    field = [field, field(1)];
end
