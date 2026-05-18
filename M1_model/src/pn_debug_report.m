function report = pn_debug_report(varargin)
%PN_DEBUG_REPORT Targeted diagnostics for the slab PN plane-source setup.

    cfg = struct( ...
        'orders', [3, 5, 7], ...
        'plane_source_method', -2, ...
        'run_plane_source', true, ...
        'verbose', true);
    cfg = apply_overrides(cfg, varargin{:});
    cfg.orders = reshape(cfg.orders, 1, []);

    report = struct();
    report.legendre = legendre_checks(cfg.orders);
    report.transport = transport_checks(cfg.orders);
    report.scattering = scattering_check(max(cfg.orders));
    report.time_integration = time_integration_check();

    case_cfg = plane_source_master_case();
    report.dirac = dirac_check(case_cfg);
    report.reference_resolution = reference_resolution_check(case_cfg);

    if cfg.run_plane_source
        ref = sn_reference_plane_source(case_cfg.reference);
        plane_source = cell(1, numel(cfg.orders));
        for idx = 1:numel(cfg.orders)
            plane_source{idx} = plane_source_check(cfg.orders(idx), cfg.plane_source_method, case_cfg, ref);
        end
        report.plane_source = [plane_source{:}];
    else
        report.plane_source = struct([]);
    end

    if cfg.verbose
        print_report(report, cfg);
    end
end

function out = legendre_checks(orders)
    out = cell(1, numel(orders));
    for idx = 1:numel(orders)
        N = orders(idx);
        [xi, w, R, L] = pn_basis(N);
        P = legpoly_eval(xi, N);
        moment_scale = (2 .* (0:N) + 1) / 2;
        expected_L = P.' * diag(moment_scale);
        expected_R = P * diag(w);

        probe = randn(N + 1, 4);
        psi_from_L = L * probe;
        psi_from_formula = expected_L * probe;

        item = struct();
        item.order = N;
        item.L_error = norm(L - expected_L, inf);
        item.R_error = norm(R - expected_R, inf);
        item.reconstruction_error = norm(psi_from_L - psi_from_formula, inf);
        item.rho_recovery_error = max(abs((w.' * psi_from_L) - probe(1, :)));
        item.passed = item.L_error < 1.0e-12 ...
            && item.R_error < 1.0e-12 ...
            && item.reconstruction_error < 1.0e-12 ...
            && item.rho_recovery_error < 1.0e-12;
        out{idx} = item;
    end
    out = [out{:}];
end

function out = transport_checks(orders)
    out = cell(1, numel(orders));
    for idx = 1:numel(orders)
        N = orders(idx);
        [xi, ~, R, L] = pn_basis(N);
        A = R * diag(xi) * L;
        A_expected = expected_transport_matrix(N);
        eigvals = sort(real(eig(A)));

        item = struct();
        item.order = N;
        item.matrix_error = norm(A - A_expected, inf);
        item.eigenvalues = eigvals(:).';
        item.passed = item.matrix_error < 1.0e-12;
        out{idx} = item;
    end
    out = [out{:}];
end

function out = scattering_check(N)
    dt = 0.137;
    sigma_a = zeros(1, 4);
    sigma_s = [0.1, 0.7, 1.3, 2.4];
    q = zeros(1, 4);
    u_in = randn(N + 1, 4);
    u_out = pn_relax_isotropic(u_in, N, dt, sigma_a, sigma_s, q);

    expected = u_in;
    expected(2:end, :) = u_in(2:end, :) ./ (ones(N, 1) * (1 + dt * sigma_s));

    out = struct();
    out.order = N;
    out.u0_error = max(abs(u_out(1, :) - u_in(1, :)));
    out.higher_moment_error = max(abs(u_out(2:end, :) - expected(2:end, :)), [], 'all');
    out.passed = out.u0_error < 1.0e-12 && out.higher_moment_error < 1.0e-12;
end

function out = time_integration_check()
    N = 7;
    dt = 0.01;
    dz = 1.0;
    sigma_a = 0.0;
    sigma_s = 1.0;
    q = 0.0;

    u = zeros(N + 1, 1);
    u(2) = 1.0;
    collision_only = pn_heun_step(u, dt, dz, N, -2, sigma_a, sigma_s, q, 'periodic', 0.0);

    exact_factor = exp(-sigma_s * dt);
    implicit_euler_factor = 1.0 / (1.0 + sigma_s * dt);
    current_factor = collision_only(2);

    case_cfg = plane_source_master_case('sigma_s_value', 0.0);
    N_transport = 31;
    u_heun = zeros(N_transport + 1, numel(case_cfg.rho0));
    u_heun(1, :) = case_cfg.rho0;
    u_single_step = u_heun;

    for step = 1:case_cfg.num_steps
        flux = pn_interface_flux(u_single_step, N_transport, -2, case_cfg.dt, case_cfg.dz, ...
            'vacuum', case_cfg.psi_vac);
        rhs = -(flux(:, 2:end) - flux(:, 1:end-1)) / case_cfg.dz;
        u_single_step = u_single_step + case_cfg.dt * rhs;

        u_heun = pn_heun_step(u_heun, case_cfg.dt, case_cfg.dz, N_transport, -2, ...
            case_cfg.sigma_a, case_cfg.sigma_s, case_cfg.source_strength, ...
            'vacuum', case_cfg.psi_vac);
    end

    obs_heun = pn_observables(u_heun);
    obs_single_step = pn_observables(u_single_step);

    out = struct();
    out.collision_dt = dt;
    out.collision_sigma_s = sigma_s;
    out.current_collision_factor = current_factor;
    out.implicit_euler_factor = implicit_euler_factor;
    out.exact_collision_factor = exact_factor;
    out.collision_factor_error_vs_ie = current_factor - implicit_euler_factor;
    out.collision_factor_error_vs_exact = current_factor - exact_factor;
    out.transport_order = N_transport;
    out.transport_linf_difference = max(abs(obs_heun.rho - obs_single_step.rho));
    out.transport_peak_heun = max(obs_heun.rho);
    out.transport_peak_single_step = max(obs_single_step.rho);
end

function out = dirac_check(case_cfg)
    baseline = 2 * case_cfg.psi_vac;
    increment = case_cfg.rho0 - baseline;
    active_cells = find(abs(increment) > 1.0e-12);

    out = struct();
    out.num_active_cells = numel(active_cells);
    out.active_cells = active_cells(:).';
    out.increment_values = increment(active_cells).';
    out.dirac_mass = sum(increment) * case_cfg.dz;
    out.total_mass = sum(case_cfg.rho0) * case_cfg.dz;
    out.passed = out.num_active_cells == 2 && abs(out.dirac_mass - 2.0) < 1.0e-12;
end

function out = reference_resolution_check(case_cfg)
    coarse_ref = sn_reference_plane_source(case_cfg.reference);
    coarse_interp = interp1(coarse_ref.z, coarse_ref.rho, case_cfg.z, 'linear', 'extrap');

    fine_reference_cfg = case_cfg.reference;
    fine_reference_cfg.n_cells = 2400;
    fine_reference_cfg.n_mu = 256;
    fine_ref = sn_reference_plane_source(fine_reference_cfg);
    fine_interp = interp1(fine_ref.z, fine_ref.rho, case_cfg.z, 'linear', 'extrap');

    out = struct();
    out.coarse_native_max = max(coarse_ref.rho);
    out.coarse_interp_max = max(coarse_interp);
    out.fine_native_max = max(fine_ref.rho);
    out.fine_interp_max = max(fine_interp);
    out.interp_linf_difference = max(abs(coarse_interp - fine_interp));
    out.coarse_interp_mass = sum(coarse_interp) * case_cfg.dz;
    out.fine_interp_mass = sum(fine_interp) * case_cfg.dz;
end

function out = plane_source_check(N, method, case_cfg, ref)
    uPN = simulate_pn(case_cfg.rho0, case_cfg.dt, case_cfg.dz, case_cfg.num_steps, ...
        N, method, case_cfg.sigma_a, case_cfg.sigma_s, case_cfg.source_strength, case_cfg.vacuum_opts);
    obs = pn_observables(uPN);
    rho_ref = interp1(ref.z, ref.rho, case_cfg.z, 'linear', 'extrap');

    speeds = positive_characteristic_speeds(N);
    wave_mass = wave_masses(case_cfg.z, obs.rho, rho_ref, speeds, case_cfg.T, case_cfg.dz);

    out = struct();
    out.order = N;
    out.method = method;
    out.total_mass = sum(obs.rho) * case_cfg.dz;
    out.reference_total_mass = sum(rho_ref) * case_cfg.dz;
    out.max_rho = max(obs.rho);
    out.reference_max_rho = max(rho_ref);
    out.wave_mass = wave_mass;
end

function wave_mass = wave_masses(z, rho, rho_ref, speeds, tf, dz)
    centers = speeds(:) * tf;
    wave_mass = cell(1, numel(centers));
    for idx = 1:numel(centers)
        center = centers(idx);
        half_width = characteristic_window(centers, idx, dz);
        mask = abs(abs(z(:)) - center) <= half_width;

        item = struct();
        item.speed = speeds(idx);
        item.center = center;
        item.half_width = half_width;
        item.pn_mass = sum(rho(mask)) * dz;
        item.reference_mass = sum(rho_ref(mask)) * dz;
        item.pn_peak = max(rho(mask));
        item.reference_peak = max(rho_ref(mask));
        wave_mass{idx} = item;
    end
    wave_mass = [wave_mass{:}];
end

function half_width = characteristic_window(centers, idx, dz)
    if numel(centers) == 1
        nearest_distance = inf;
    elseif idx == 1
        nearest_distance = centers(2) - centers(1);
    elseif idx == numel(centers)
        nearest_distance = centers(end) - centers(end - 1);
    else
        nearest_distance = min(centers(idx + 1) - centers(idx), centers(idx) - centers(idx - 1));
    end

    if isfinite(nearest_distance)
        half_width = max(2 * dz, 0.45 * nearest_distance);
    else
        half_width = 2 * dz;
    end
end

function A = expected_transport_matrix(N)
    A = zeros(N + 1, N + 1);
    for ell = 0:N
        row = ell + 1;
        if ell < N
            A(row, row + 1) = (ell + 1) / (2 * ell + 1);
        end
        if ell > 0
            A(row, row - 1) = ell / (2 * ell + 1);
        end
    end
end

function speeds = positive_characteristic_speeds(N)
    A = expected_transport_matrix(N);
    eigvals = sort(real(eig(A)));
    speeds = eigvals(eigvals > 1.0e-12).';
end

function print_report(report, cfg)
    fprintf('PN debug report (orders = %s, method = %s)\n', ...
        mat2str(cfg.orders), method_name(cfg.plane_source_method));

    fprintf('\n1. Legendre normalization\n');
    for idx = 1:numel(report.legendre)
        item = report.legendre(idx);
        fprintf(['  P%d: %s, L err = %.3e, R err = %.3e, ' ...
            'reconstruction err = %.3e, rho recovery err = %.3e\n'], ...
            item.order, status_text(item.passed), item.L_error, item.R_error, ...
            item.reconstruction_error, item.rho_recovery_error);
    end

    fprintf('\n2. PN transport matrix\n');
    for idx = 1:numel(report.transport)
        item = report.transport(idx);
        fprintf('  P%d: %s, matrix err = %.3e, eigenvalues = %s\n', ...
            item.order, status_text(item.passed), item.matrix_error, mat2str(item.eigenvalues, 8));
    end

    fprintf('\n3. Isotropic scattering\n');
    fprintf('  P%d: %s, u0 err = %.3e, higher-moment err = %.3e\n', ...
        report.scattering.order, status_text(report.scattering.passed), ...
        report.scattering.u0_error, report.scattering.higher_moment_error);

    fprintf('\n4. Time integration\n');
    fprintf(['  collision-only factor at dt = %.3g, sigma_s = %.3g: current = %.12f, ' ...
        'implicit Euler = %.12f, exact = %.12f\n'], ...
        report.time_integration.collision_dt, report.time_integration.collision_sigma_s, ...
        report.time_integration.current_collision_factor, ...
        report.time_integration.implicit_euler_factor, ...
        report.time_integration.exact_collision_factor);
    fprintf(['  pure-transport P%d, LW single-step vs Heun(LW-flux): ' ...
        'rho Linf difference = %.6f, peak(single-step) = %.6f, peak(Heun) = %.6f\n'], ...
        report.time_integration.transport_order, ...
        report.time_integration.transport_linf_difference, ...
        report.time_integration.transport_peak_single_step, ...
        report.time_integration.transport_peak_heun);

    fprintf('\n5. Dirac initialization\n');
    fprintf(['  %s, active cells = %d, dirac mass = %.16g, total mass = %.16g, ' ...
        'increments = %s\n'], ...
        status_text(report.dirac.passed), report.dirac.num_active_cells, ...
        report.dirac.dirac_mass, report.dirac.total_mass, ...
        mat2str(report.dirac.increment_values, 8));

    fprintf('\n6. Reference resolution\n');
    fprintf(['  coarse SN (300x64): native max = %.8f, interp max = %.8f, interp mass = %.8f\n'], ...
        report.reference_resolution.coarse_native_max, ...
        report.reference_resolution.coarse_interp_max, ...
        report.reference_resolution.coarse_interp_mass);
    fprintf(['  fine   SN (2400x256): native max = %.8f, interp max = %.8f, interp mass = %.8f\n'], ...
        report.reference_resolution.fine_native_max, ...
        report.reference_resolution.fine_interp_max, ...
        report.reference_resolution.fine_interp_mass);
    fprintf('  coarse-vs-fine interp rho Linf difference = %.8f\n', ...
        report.reference_resolution.interp_linf_difference);

    if isempty(report.plane_source)
        return;
    end

    fprintf('\n7. Peak diagnostics\n');
    for idx = 1:numel(report.plane_source)
        item = report.plane_source(idx);
        fprintf(['  P%d-%s: total mass = %.8f (ref %.8f), max rho = %.8f (ref %.8f)\n'], ...
            item.order, method_name(item.method), item.total_mass, item.reference_total_mass, ...
            item.max_rho, item.reference_max_rho);
        for wave_idx = 1:numel(item.wave_mass)
            wave = item.wave_mass(wave_idx);
            fprintf(['    lambda = %.6f: mass(|z|-lambda t <= %.4f) = %.8f (ref %.8f), ' ...
                'peak = %.8f (ref %.8f)\n'], ...
                wave.speed, wave.half_width, wave.pn_mass, wave.reference_mass, ...
                wave.pn_peak, wave.reference_peak);
        end
    end
end

function name = method_name(method)
    switch method
        case 1
            name = 'MCL';
        case -1
            name = 'LLF';
        case -2
            name = 'LW';
        otherwise
            name = sprintf('method=%g', method);
    end
end

function text = status_text(passed)
    if passed
        text = 'OK';
    else
        text = 'VIOLATION';
    end
end

function cfg = apply_overrides(cfg, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('pn_debug_report:invalidOverrides', ...
            'Overrides must be provided as name/value pairs.');
    end

    for idx = 1:2:numel(varargin)
        name = varargin{idx};
        value = varargin{idx + 1};
        if ~isfield(cfg, name)
            error('pn_debug_report:unknownField', ...
                'Unknown override field "%s".', name);
        end
        cfg.(name) = value;
    end
end
