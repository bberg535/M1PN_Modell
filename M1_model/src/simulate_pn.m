function uPN = simulate_pn(rho0, dt, dz, num_steps, N, method, sigma_a, sigma_s, q, opts)
%SIMULATE_PN Advance the PN model with SSP Heun time stepping.
% Boundary handling is shared with the active HOLO coupling so that pure PN
% and coupled PN use the same transport conventions.

    if nargin < 10 || ~isstruct(opts)
        opts = struct();
    end

    boundary = get_field_or(opts, 'boundary', 'periodic');
    psi_boundary = get_field_or(opts, 'psi_boundary', 0.0);
    nec = numel(rho0);

    sigma_a = reshape_cell_field(sigma_a, nec, 'simulate_pn:sigmaA');
    sigma_s = reshape_cell_field(sigma_s, nec, 'simulate_pn:sigmaS');
    q = reshape_cell_field(q, nec, 'simulate_pn:source');

    uPN = zeros(N + 1, nec);
    uPN(1, :) = reshape(rho0, 1, []);

    for step = 1:num_steps
        k1 = transport_rhs(uPN, N, method, dt, dz, boundary, psi_boundary);
        uPN1 = uPN + dt * k1;
        uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, q);

        k2 = transport_rhs(uPN1, N, method, dt, dz, boundary, psi_boundary);
        uPN2 = uPN1 + dt * k2;
        uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, q);

        uPN = 0.5 * (uPN + uPN2);

        if any(~isfinite(uPN(:))) || any(~isreal(uPN(:)))
            error('simulate_pn:invalidState', ...
                'PN state became invalid at step %d.', step);
        end
    end
end

function rhs = transport_rhs(uPN, N, method, dt, dz, boundary, psi_boundary)
    flux = pn_interface_flux(uPN, N, method, dt, dz, boundary, psi_boundary);
    rhs = -(flux(:, 2:end) - flux(:, 1:end-1)) / dz;
end

function field = reshape_cell_field(field, nec, error_id)
    if isscalar(field)
        field = field * ones(1, nec);
    else
        field = reshape(field, 1, []);
    end

    if numel(field) ~= nec
        error(error_id, 'Expected field length %d.', nec);
    end
end

function value = get_field_or(opts, name, default_value)
    if isfield(opts, name)
        value = opts.(name);
    else
        value = default_value;
    end
end
