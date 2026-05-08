function out = simulate_pn_siu(cfg)
%SIMULATE_PN_SIU Semi-implicit upwind PN solver from Positive PN Closures.
% Implements the slab-geometry SIU scheme from Hauck and McClarren,
% "Positive PN Closures", SIAM J. Sci. Comput. 32(5), 2010:
% - transport-relaxation form (13),
% - SIU update (15)-(20),
% - periodic plane-source example from Section 3.2 / Figure 3.

    if nargin < 1
        cfg = struct();
    end

    N = get_field_or(cfg, 'N', 3);
    dt = get_field_or(cfg, 'dt', 1.0e-3);
    dz = get_field_or(cfg, 'dz', 2.0e-3);
    tf = get_field_or(cfg, 'tf', 1.0);
    domain_length = get_field_or(cfg, 'domain_length', 2.0);
    sigma = get_field_or(cfg, 'sigma', 1.0);
    mu_eval = get_field_or(cfg, 'mu_eval', linspace(-1.0, 1.0, 1000));

    num_cells = round(domain_length / dz);
    if abs(num_cells * dz - domain_length) > 1e-12 * max(1, domain_length)
        error('simulate_pn_siu:meshMismatch', ...
            'Domain length %.16g is not an integer multiple of dz=%.16g.', ...
            domain_length, dz);
    end

    num_steps = round(tf / dt);
    if abs(num_steps * dt - tf) > 1e-12 * max(1, tf)
        error('simulate_pn_siu:timeMismatch', ...
            'tf=%.16g is not an integer multiple of dt=%.16g.', tf, dt);
    end

    [lambda, w, R, L] = pn_basis(N);
    G = diag(2 ./ (2 .* (0:N) + 1));

    Lambda = diag(lambda);
    lambda_plus = max(Lambda, 0);
    lambda_minus = min(Lambda, 0);
    C = eye(N + 1) - (dt / dz) * (lambda_plus - lambda_minus);
    Cplus = -(dt / dz) * lambda_minus;
    Cminus = (dt / dz) * lambda_plus;

    gamma = 1.0 / (1.0 + sigma * dt);

    u = zeros(N + 1, num_cells);
    u(1, 1) = 1.0 / dz;
    v = L * u;

    left_shift = [num_cells, 1:(num_cells - 1)];
    right_shift = [2:num_cells, 1];

    for step = 1:num_steps
        cd_state = C * v + Cminus * v(:, left_shift) + Cplus * v(:, right_shift);
        rho = w.' * cd_state;
        v = gamma * cd_state + (1.0 - gamma) * (rho / 2.0);

        if any(~isfinite(v(:))) || any(~isreal(v(:)))
            error('simulate_pn_siu:invalidState', ...
                'SIU state became invalid at step %d.', step);
        end
    end

    u = R * v;
    coeff = G \ u;
    p_eval = legpoly_eval(mu_eval, N);
    reconstruction = coeff.' * p_eval;
    min_reconstruction = min(reconstruction, [], 2);
    [~, worst_idx] = min(min_reconstruction);

    z = (dz / 2.0):dz:(domain_length - dz / 2.0);
    rho = w.' * v;

    out = struct();
    out.N = N;
    out.dt = dt;
    out.dz = dz;
    out.tf = tf;
    out.sigma = sigma;
    out.z = z(:);
    out.u = u;
    out.v = v;
    out.rho = rho(:);
    out.lambda = lambda(:);
    out.w = w(:);
    out.mu_eval = mu_eval(:);
    out.reconstruction = reconstruction;
    out.min_reconstruction = min_reconstruction(:);
    out.worst_index = worst_idx;
    out.worst_z = z(worst_idx);
    out.worst_profile = reconstruction(worst_idx, :).';
    out.quadrature_values = v(:, worst_idx);
    out.max_ratio = max(abs(min_reconstruction) ./ max(rho(:) / 2.0, 1.0e-14));
end

function value = get_field_or(s, name, default_value)
    if isfield(s, name)
        value = s.(name);
    else
        value = default_value;
    end
end
