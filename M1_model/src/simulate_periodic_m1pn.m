function [u, uPN] = simulate_periodic_m1pn(rho0, dt, dz, num_steps, NPN, ilim, sigma_a, sigma_s, source_strength)
%SIMULATE_PERIODIC_M1PN Compatibility wrapper for periodic coupled runs.

    if ilim ~= -1
        error('simulate_periodic_m1pn:unsupportedLegacySwitch', ...
            ['The legacy periodic M1PN wrapper only accepts ilim = -1. ' ...
            'Configure the PN HO candidate via simulate_m1pn(..., opts).']);
    end

    opts = struct('boundary', 'periodic', 'psi_boundary', 0.0, 'pn_ho_method', 1);
    [u_cell, uPN_cell] = simulate_m1pn(rho0, dt, dz, num_steps, NPN, ...
        sigma_a, sigma_s, source_strength, opts);
    u = [u_cell, u_cell(:, 1)];
    uPN = [uPN_cell, uPN_cell(:, 1)];
end
