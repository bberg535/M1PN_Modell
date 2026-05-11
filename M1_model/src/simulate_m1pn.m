function [u, uPN, diagnostics] = simulate_m1pn(rho0, dt, dz, num_steps, pn_order, sigma_a, sigma_s, source_strength, opts)
%SIMULATE_M1PN Advance the coupled HOLO M1PN model.
% The returned M1 state is the actual model solution. The PN state is an
% auxiliary high-order transport state used for the HO candidate flux.
%
% Available coupling modes via opts.coupling_mode:
% - explicit_failsafe (default): current explicit HO correction path
% - explicit_synced: current path plus least-squares PN moment sync
% - implicit_holo: outer Picard HOLO iteration with safe fallback

    if nargin < 9 || ~isstruct(opts)
        opts = struct();
    end
    solver_cfg = m1pn_solver_config(opts);

    nec = numel(rho0);
    sigma_a = reshape_cell_field(sigma_a, nec, 'simulate_m1pn:sigmaA');
    sigma_s = reshape_cell_field(sigma_s, nec, 'simulate_m1pn:sigmaS');

    u = zeros(2, nec);
    uPN = zeros(pn_order + 1, nec);
    u(1, :) = reshape(rho0, 1, []);
    uPN(1, :) = reshape(rho0, 1, []);

    diagnostics = struct();
    diagnostics.coupling_mode = solver_cfg.coupling_mode;
    diagnostics.time = zeros(num_steps, 1);
    diagnostics.min_alpha = ones(num_steps, 1);
    diagnostics.max_alpha = ones(num_steps, 1);
    diagnostics.num_limited = zeros(num_steps, 1);
    diagnostics.outer_iterations = zeros(num_steps, 1);
    diagnostics.inner_iterations = zeros(num_steps, 1);
    diagnostics.fallback_used = false(num_steps, 1);
    diagnostics.consistency_norm = zeros(num_steps, 1);
    diagnostics.sync_correction = zeros(num_steps, 1);

    for step = 1:num_steps
        step_info = struct('mode', 'M1PN', 'step', step, 'time', step * dt);
        [u, uPN, step_diag] = advance_step(u, uPN, dt, dz, sigma_a, sigma_s, ...
            pn_order, source_strength, solver_cfg, step_info);
        diagnostics.time(step) = step * dt;
        diagnostics.min_alpha(step) = step_diag.min_alpha;
        diagnostics.max_alpha(step) = step_diag.max_alpha;
        diagnostics.num_limited(step) = step_diag.num_limited;
        diagnostics.outer_iterations(step) = get_diag_field(step_diag, 'outer_iterations', 0);
        diagnostics.inner_iterations(step) = get_diag_field(step_diag, 'inner_iterations', 0);
        diagnostics.fallback_used(step) = get_diag_field(step_diag, 'fallback_used', false);
        diagnostics.consistency_norm(step) = get_diag_field(step_diag, 'consistency_norm', 0.0);
        diagnostics.sync_correction(step) = get_diag_field(step_diag, 'max_sync_correction', 0.0);
    end
end

function [u, uPN, step_diag] = advance_step(u, uPN, dt, dz, sigma_a, sigma_s, ...
        pn_order, source_strength, solver_cfg, step_info)
    switch solver_cfg.coupling_mode
        case 'explicit_failsafe'
            [u, uPN, step_diag] = heuns_method(u, uPN, dt, dz, sigma_a, sigma_s, ...
                pn_order, source_strength, solver_cfg, step_info);
            if solver_cfg.sync_auxiliary_pn
                [uPN, sync_diag] = pn_sync_low_moments(uPN, u(1:2, :));
                step_diag = add_sync_diag(step_diag, sync_diag);
            else
                step_diag = add_sync_diag(step_diag, []);
            end

        case 'explicit_synced'
            [u, uPN, step_diag] = heuns_method(u, uPN, dt, dz, sigma_a, sigma_s, ...
                pn_order, source_strength, solver_cfg, step_info);
            [uPN, sync_diag] = pn_sync_low_moments(uPN, u(1:2, :));
            step_diag = add_sync_diag(step_diag, sync_diag);

        case 'implicit_holo'
            [u, uPN, step_diag] = m1pn_implicit_holo_step(u, uPN, dt, dz, ...
                sigma_a, sigma_s, pn_order, source_strength, solver_cfg, step_info);

        otherwise
            error('simulate_m1pn:unsupportedCouplingMode', ...
                'Unsupported coupling mode "%s".', solver_cfg.coupling_mode);
    end
end

function step_diag = add_sync_diag(step_diag, sync_diag)
    if isempty(sync_diag)
        step_diag.max_sync_correction = 0.0;
        step_diag.max_pre_sync_mismatch = 0.0;
        step_diag.max_post_sync_mismatch = 0.0;
        return;
    end

    step_diag.max_sync_correction = sync_diag.max_correction;
    step_diag.max_pre_sync_mismatch = sync_diag.max_pre_mismatch;
    step_diag.max_post_sync_mismatch = sync_diag.max_post_mismatch;
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

function value = get_diag_field(diag_struct, field_name, default_value)
    if isfield(diag_struct, field_name)
        value = diag_struct.(field_name);
    else
        value = default_value;
    end
end
