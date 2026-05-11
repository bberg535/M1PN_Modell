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
        uPN = pn_heun_step(uPN, dt, dz, N, method, sigma_a, sigma_s, q, ...
            boundary, psi_boundary);

        if any(~isfinite(uPN(:))) || any(~isreal(uPN(:)))
            error('simulate_pn:invalidState', ...
                'PN state became invalid at step %d.', step);
        end
    end
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
