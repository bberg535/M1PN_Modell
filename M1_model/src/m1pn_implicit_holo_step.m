function [u, uPN, step_diag] = m1pn_implicit_holo_step(u_old, uPN_old, dt, dz, ...
        sigma_a, sigma_s, N, source_strength, solver_cfg, step_info)
%M1PN_IMPLICIT_HOLO_STEP Picard-style implicit HOLO step with safe fallback.
% The nonlinear iteration is driven by the LO M1 state. The PN state is
% nonlinearly eliminated by synchronizing a raw HO predictor onto the LO
% moments and by feeding the resulting consistency flux back into the LO
% solve.

    if size(u_old, 2) ~= size(uPN_old, 2)
        error('m1pn_implicit_holo_step:dimensionMismatch', ...
            'Expected M1 and PN states to use the same number of spatial cells.');
    end

    boundary = solver_cfg.boundary;
    psi_boundary = solver_cfg.psi_boundary;
    pn_ho_method = solver_cfg.pn_ho_method;
    source = source_term(size(u_old, 2), source_strength);
    pn_source = source(1, :);

    predictor_opts = struct('boundary', boundary, ...
        'psi_boundary', psi_boundary, ...
        'pn_ho_method', pn_ho_method);
    [u_predict, uPN_predict, predictor_diag] = heuns_method(u_old, uPN_old, ...
        dt, dz, sigma_a, sigma_s, N, source_strength, predictor_opts, step_info);
    [uPN_predict, predictor_sync_diag] = pn_sync_low_moments(uPN_predict, u_predict(1:2, :));
    uPN_raw = pn_heun_step(uPN_old, dt, dz, N, pn_ho_method, sigma_a, ...
        sigma_s, pn_source, boundary, psi_boundary);
    if any(~isfinite(uPN_raw(:))) || any(~isreal(uPN_raw(:)))
        u = u_predict;
        uPN = uPN_predict;
        step_diag = build_fallback_diag(predictor_diag, 0, 0, inf, inf, ...
            predictor_sync_diag);
        assert_valid_m1pn_state(u, uPN, step_info, step_diag);
        return;
    end

    switch solver_cfg.holo_initial_guess
        case 'previous_state'
            u_guess = u_old;
        otherwise
            u_guess = u_predict;
    end

    outer_converged = false;
    last_rel_change = inf;
    last_consistency_norm = inf;
    max_inner_iterations = 0;
    final_sync_diag = zero_sync_diag();
    outer_iter = 0;

    for outer_iter = 1:solver_cfg.holo_outer_max_iters
        [uPN_sync, sync_diag] = pn_sync_low_moments(uPN_raw, u_guess(1:2, :));
        flux_lo_lag = m1_interface_flux(u_guess, -1, dt, dz, boundary, psi_boundary);
        flux_ho = pn_interface_flux(uPN_sync, N, pn_ho_method, dt, dz, ...
            boundary, psi_boundary);
        consistency_flux = flux_ho(1:2, :) - flux_lo_lag;
        last_consistency_norm = max(abs(consistency_flux(:)));

        [u_candidate, inner_diag] = solve_lo_state(u_old, u_guess, ...
            consistency_flux, dt, dz, source, sigma_a, sigma_s, ...
            boundary, psi_boundary, solver_cfg);
        max_inner_iterations = max(max_inner_iterations, inner_diag.iterations);

        if ~inner_diag.valid_state
            final_sync_diag = sync_diag;
            break;
        end

        u_candidate = solver_cfg.holo_relaxation * u_candidate + ...
            (1.0 - solver_cfg.holo_relaxation) * u_guess;
        last_rel_change = relative_change(u_candidate, u_guess);
        u_guess = u_candidate;
        final_sync_diag = sync_diag;

        if last_rel_change <= solver_cfg.holo_tol
            outer_converged = true;
            break;
        end
    end

    [uPN_final, final_sync_diag] = pn_sync_low_moments(uPN_raw, u_guess(1:2, :));
    final_state_valid = is_valid_coupled_state(u_guess, uPN_final);
    if outer_converged && final_state_valid
        u = u_guess;
        uPN = uPN_final;
        step_diag = build_implicit_diag(size(u_old, 2) + 1, outer_iter, ...
            max_inner_iterations, false, last_rel_change, ...
            last_consistency_norm, final_sync_diag);
    else
        u = u_predict;
        uPN = uPN_predict;
        step_diag = build_fallback_diag(predictor_diag, outer_iter, ...
            max_inner_iterations, last_rel_change, last_consistency_norm, ...
            combine_sync_diag(final_sync_diag, predictor_sync_diag));
    end

    assert_valid_m1pn_state(u, uPN, step_info, step_diag);
end

function [u_next, diag] = solve_lo_state(u_old, u_init, consistency_flux, dt, dz, ...
        source, sigma_a, sigma_s, boundary, psi_boundary, solver_cfg)
    u_iter = u_init;
    diag = struct('iterations', 0, 'converged', false, ...
        'last_rel_change', inf, 'valid_state', true);

    for iter = 1:solver_cfg.holo_inner_max_iters
        if ~is_valid_lo_state(u_iter)
            diag.iterations = iter - 1;
            diag.valid_state = false;
            u_next = u_iter;
            return;
        end

        flux_lo = m1_interface_flux(u_iter, -1, dt, dz, boundary, psi_boundary);
        flux_total = flux_lo + consistency_flux;
        transport = -(flux_total(:, 2:end) - flux_total(:, 1:end-1)) / dz;

        rhs = u_old + dt * (transport + source);
        u_next = zeros(size(u_iter));
        u_next(1, :) = rhs(1, :) ./ (1 + dt * sigma_a);
        u_next(2, :) = rhs(2, :) ./ (1 + dt * (sigma_a + sigma_s));
        u_next = enforce_lo_realizability(u_iter, u_next);

        diag.iterations = iter;
        diag.last_rel_change = relative_change(u_next, u_iter);
        u_iter = u_next;

        if diag.last_rel_change <= solver_cfg.holo_tol
            diag.converged = true;
            break;
        end
    end
end

function step_diag = build_implicit_diag(interface_count, outer_iterations, ...
        inner_iterations, fallback_used, last_rel_change, ...
        consistency_norm, sync_diag)
    step_diag = struct();
    step_diag.alpha = ones(1, interface_count);
    step_diag.min_alpha = 1.0;
    step_diag.max_alpha = 1.0;
    step_diag.num_limited = 0;
    step_diag.outer_iterations = outer_iterations;
    step_diag.inner_iterations = inner_iterations;
    step_diag.fallback_used = fallback_used;
    step_diag.last_rel_change = last_rel_change;
    step_diag.consistency_norm = consistency_norm;
    step_diag.max_sync_correction = sync_diag.max_correction;
    step_diag.max_pre_sync_mismatch = sync_diag.max_pre_mismatch;
    step_diag.max_post_sync_mismatch = sync_diag.max_post_mismatch;
end

function step_diag = build_fallback_diag(predictor_diag, outer_iterations, ...
        inner_iterations, last_rel_change, consistency_norm, sync_diag)
    step_diag = predictor_diag;
    step_diag.outer_iterations = outer_iterations;
    step_diag.inner_iterations = inner_iterations;
    step_diag.fallback_used = true;
    step_diag.last_rel_change = last_rel_change;
    step_diag.consistency_norm = consistency_norm;
    step_diag.max_sync_correction = sync_diag.max_correction;
    step_diag.max_pre_sync_mismatch = sync_diag.max_pre_mismatch;
    step_diag.max_post_sync_mismatch = sync_diag.max_post_mismatch;
end

function ok = is_valid_lo_state(u)
    if any(~isfinite(u(:))) || any(~isreal(u(:)))
        ok = false;
        return;
    end

    rho = u(1, :);
    j = u(2, :);
    tol = 1.0e-10;
    ok = all(rho >= -tol) && all(abs(j) <= rho + tol);
end

function ok = is_valid_coupled_state(u, uPN)
    ok = is_valid_lo_state(u) && all(isfinite(uPN(:))) && all(isreal(uPN(:)));
end

function u_limited = enforce_lo_realizability(u_base, u_trial)
    rho0 = u_base(1, :);
    j0 = u_base(2, :);
    rho1 = u_trial(1, :);
    j1 = u_trial(2, :);

    q0 = [rho0; rho0 + j0; rho0 - j0];
    q1 = [rho1; rho1 + j1; rho1 - j1];

    theta = ones(1, size(u_base, 2));
    tol = 1.0e-12;
    for row = 1:size(q0, 1)
        invalid_target = q1(row, :) < tol;
        needs_limit = invalid_target & (q0(row, :) > tol);
        theta(needs_limit) = min(theta(needs_limit), ...
            q0(row, needs_limit) ./ max(q0(row, needs_limit) - q1(row, needs_limit), tol));
        theta(invalid_target & ~(q0(row, :) > tol)) = 0.0;
    end

    theta = max(0.0, min(1.0, 0.999 * theta));
    u_limited = u_base + (u_trial - u_base) .* theta(ones(size(u_base, 1), 1), :);

    rho = max(u_limited(1, :), 0.0);
    j = max(min(u_limited(2, :), rho), -rho);
    u_limited = [rho; j];
end

function value = relative_change(new_state, old_state)
    scale = max(1.0, max(abs(new_state(:))));
    value = max(abs(new_state(:) - old_state(:))) / scale;
end

function sync_diag = zero_sync_diag()
    sync_diag = struct();
    sync_diag.max_pre_mismatch = 0.0;
    sync_diag.max_post_mismatch = 0.0;
    sync_diag.max_correction = 0.0;
end

function sync_diag = combine_sync_diag(primary_diag, secondary_diag)
    sync_diag = struct();
    sync_diag.max_pre_mismatch = max(primary_diag.max_pre_mismatch, secondary_diag.max_pre_mismatch);
    sync_diag.max_post_mismatch = max(primary_diag.max_post_mismatch, secondary_diag.max_post_mismatch);
    sync_diag.max_correction = max(primary_diag.max_correction, secondary_diag.max_correction);
end
