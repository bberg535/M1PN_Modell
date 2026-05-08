function u = simulate_m1(rho0, dt, dz, num_steps, method, sigma_a, sigma_s, source_strength, opts)
%SIMULATE_M1 Advance the closed M1 model with SSP Heun time stepping.
% This is the reference low-order/standalone M1 transport solver used by
% the HOLO coupling and by the standalone benchmarks.

    if nargin < 9 || ~isstruct(opts)
        opts = struct();
    end

    boundary = get_field_or(opts, 'boundary', 'periodic');
    psi_boundary = get_field_or(opts, 'psi_boundary', 0.0);
    nec = numel(rho0);

    sigma_a = reshape_cell_field(sigma_a, nec, 'simulate_m1:sigmaA');
    sigma_s = reshape_cell_field(sigma_s, nec, 'simulate_m1:sigmaS');
    source = source_term(nec, source_strength);

    u = zeros(2, nec);
    u(1, :) = reshape(rho0, 1, []);

    for step = 1:num_steps
        rhs1 = build_rhs(u, method, dt, dz, boundary, psi_boundary, ...
            source, sigma_a, sigma_s);
        u1 = u + dt * rhs1;

        rhs2 = build_rhs(u1, method, dt, dz, boundary, psi_boundary, ...
            source, sigma_a, sigma_s);
        u2 = u1 + dt * rhs2;

        u = 0.5 * (u + u2);

        if any(~isfinite(u(:))) || any(~isreal(u(:)))
            error('simulate_m1:invalidState', ...
                'M1 state became invalid at step %d.', step);
        end
    end
end

function rhs = build_rhs(u, method, dt, dz, boundary, psi_boundary, source, sigma_a, sigma_s)
    flux = m1_interface_flux(u, method, dt, dz, boundary, psi_boundary);
    rhs = -(flux(:, 2:end) - flux(:, 1:end-1)) / dz ...
        + source - reaction(u, sigma_a, sigma_s);
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
