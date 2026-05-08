function [u, uPN, diagnostics] = simulate_m1pn(rho0, dt, dz, num_steps, pn_order, sigma_a, sigma_s, source_strength, opts)
%SIMULATE_M1PN Advance the coupled HOLO M1PN model.
% The returned M1 state is the actual model solution. The PN state is an
% auxiliary high-order transport state used for the HO candidate flux.

    if nargin < 9 || ~isstruct(opts)
        opts = struct();
    end

    nec = numel(rho0);
    sigma_a = reshape_cell_field(sigma_a, nec, 'simulate_m1pn:sigmaA');
    sigma_s = reshape_cell_field(sigma_s, nec, 'simulate_m1pn:sigmaS');

    u = zeros(2, nec);
    uPN = zeros(pn_order + 1, nec);
    u(1, :) = reshape(rho0, 1, []);
    uPN(1, :) = reshape(rho0, 1, []);

    diagnostics = struct();
    diagnostics.time = zeros(num_steps, 1);
    diagnostics.min_alpha = ones(num_steps, 1);
    diagnostics.max_alpha = ones(num_steps, 1);
    diagnostics.num_limited = zeros(num_steps, 1);

    for step = 1:num_steps
        step_info = struct('mode', 'M1PN', 'step', step, 'time', step * dt);
        [u, uPN, step_diag] = heuns_method(u, uPN, dt, dz, sigma_a, sigma_s, ...
            pn_order, source_strength, opts, step_info);
        diagnostics.time(step) = step * dt;
        diagnostics.min_alpha(step) = step_diag.min_alpha;
        diagnostics.max_alpha(step) = step_diag.max_alpha;
        diagnostics.num_limited(step) = step_diag.num_limited;
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
