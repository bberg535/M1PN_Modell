function [u, uPN, step_diag] = heuns_method(u, uPN, dt, dz, sigma_a, sigma_s, N, source_strength, opts, step_info)
%HEUNS_METHOD One SSP-Heun step for the coupled HOLO M1PN model.
% M1 stays the resolved model for (rho, j). PN only supplies the
% high-order candidate transport flux used inside the fail-safe correction.

    if nargin < 9 || ~isstruct(opts)
        opts = struct();
    end
    if nargin < 10 || ~isstruct(step_info)
        step_info = struct();
    end
    if ~isfield(step_info, 'mode')
        step_info.mode = 'M1PN';
    end

    if size(u, 2) ~= size(uPN, 2)
        error('heuns_method:dimensionMismatch', ...
            'Expected u and uPN to use the same number of spatial cells.');
    end

    boundary = get_field_or(opts, 'boundary', 'periodic');
    psi_boundary = get_field_or(opts, 'psi_boundary', 0.0);
    pn_ho_method = get_field_or(opts, 'pn_ho_method', 1);
    source = source_term(size(u, 2), source_strength);
    pn_source = source(1, :);

    [rhs1_u, stage1_diag] = build_rhs(u, uPN, dt, dz, N, boundary, psi_boundary, ...
        pn_ho_method, source, sigma_a, sigma_s);
    u1 = u + dt * rhs1_u;

    pn_rhs1 = pn_transport_rhs(uPN, N, pn_ho_method, dt, dz, boundary, psi_boundary);
    uPN1 = uPN + dt * pn_rhs1;
    uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, pn_source);

    [rhs2_u, stage2_diag] = build_rhs(u1, uPN1, dt, dz, N, boundary, psi_boundary, ...
        pn_ho_method, source, sigma_a, sigma_s);
    u2 = u1 + dt * rhs2_u;

    pn_rhs2 = pn_transport_rhs(uPN1, N, pn_ho_method, dt, dz, boundary, psi_boundary);
    uPN2 = uPN1 + dt * pn_rhs2;
    uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, pn_source);

    u = 0.5 * (u + u2);
    uPN = 0.5 * (uPN + uPN2);

    step_diag = combine_stage_diagnostics(stage1_diag, stage2_diag);
    assert_valid_m1pn_state(u, uPN, step_info, step_diag);
end

function [rhs, limiter_diag] = build_rhs(u, uPN, dt, dz, N, boundary, psi_boundary, pn_ho_method, source, sigma_a, sigma_s)
    [flux_star, limiter_diag] = real_HO_flux(uPN, u, N, dt, dz, pn_ho_method, boundary, psi_boundary);
    rhs = -(flux_star(:, 2:end) - flux_star(:, 1:end-1)) / dz ...
        + source - reaction(u, sigma_a, sigma_s);
end

function rhs = pn_transport_rhs(uPN, N, method, dt, dz, boundary, psi_boundary)
    flux = pn_interface_flux(uPN, N, method, dt, dz, boundary, psi_boundary);
    rhs = -(flux(:, 2:end) - flux(:, 1:end-1)) / dz;
end

function diag_out = combine_stage_diagnostics(diag_a, diag_b)
    diag_out = struct();
    diag_out.alpha = min(diag_a.alpha, diag_b.alpha);
    diag_out.min_alpha = min(diag_a.min_alpha, diag_b.min_alpha);
    diag_out.max_alpha = max(diag_a.max_alpha, diag_b.max_alpha);
    diag_out.num_limited = max(diag_a.num_limited, diag_b.num_limited);
end

function value = get_field_or(opts, name, default_value)
    if isfield(opts, name)
        value = opts.(name);
    else
        value = default_value;
    end
end
