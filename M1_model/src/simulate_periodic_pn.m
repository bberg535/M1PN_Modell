function uPN = simulate_periodic_pn(rho0, dt, dz, num_steps, N, ilim, sigma_a, sigma_s, q)
%SIMULATE_PERIODIC_PN Compatibility wrapper for periodic PN runs.

    opts = struct('boundary', 'periodic', 'psi_boundary', 0.0);
    uPN_cell = simulate_pn(rho0, dt, dz, num_steps, N, ilim, sigma_a, sigma_s, q, opts);
    uPN = [uPN_cell, uPN_cell(:, 1)];
end
