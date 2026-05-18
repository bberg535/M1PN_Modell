% Plane-source benchmark against the S_N reference convention used in
% masterseminar_basis_vergleich for Section 6.1.1 of Seminarquelle 2.
%
% The active coupled result uses the HOLO path:
% closed low-order M1 flux + PN high-order candidate flux + fail-safe rho limiter.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

%% Benchmark setup
compare_models = struct( ...
    'mcl', false, ...
    'llf', false, ...
    'lw', false, ...
    'pn', false, ...
    'm1pn', true);

pn_ilim = -21;
pn_order = 111;
m1pn_order = 19;
% PN transport method used inside M1PN:
% - 'MCL'
% - 'LLF'
% - 'LW'
m1pn_pn_ho_method = 'MCL';
% M1PN coupling modes to compare:
% - 'explicit_failsafe' (current default path)
% - 'explicit_synced'
% - 'implicit_holo'
% Example:
% m1pn_coupling_modes = {'explicit_failsafe', 'explicit_synced', 'implicit_holo'};
m1pn_coupling_modes = {'explicit_synced'};
% Optional extra parameters for all selected M1PN coupling modes.
% Example:
% m1pn_mode_opts = struct('holo_outer_max_iters', 12, 'holo_tol', 1.0e-4);
m1pn_mode_opts = struct();

case_cfg = plane_source_master_case();
z = case_cfg.z;
rho0 = case_cfg.rho0;
dz = case_cfg.dz;
dt = case_cfg.dt;
num_cells = case_cfg.num_cells;
num_steps = case_cfg.num_steps;
sigma_a = case_cfg.sigma_a;
sigma_s = case_cfg.sigma_s;
source_strength = case_cfg.source_strength;
vacuum_opts = case_cfg.vacuum_opts;

%% Reference solution
sn_ref = sn_reference_plane_source(case_cfg.reference);

rho_ref = interp1(sn_ref.z, sn_ref.rho, z, 'linear', 'extrap');
j_ref = interp1(sn_ref.z, sn_ref.j, z, 'linear', 'extrap');
m2_ref = interp1(sn_ref.z, sn_ref.m2, z, 'linear', 'extrap');
sn_ref_moments = sn_reference_pn_moments(sn_ref, 2);
u2_ref = interp1(sn_ref.z, sn_ref_moments(3, :).', z, 'linear', 'extrap');

%% Compare models
ilim_values = [];
if compare_models.mcl
    ilim_values(end + 1) = 1;
end
if compare_models.llf
    ilim_values(end + 1) = -1;
end
if compare_models.lw
    ilim_values(end + 1) = -2;
end

m1pn_mode_list = normalize_mode_list(m1pn_coupling_modes);
num_m1pn_results = compare_models.m1pn * numel(m1pn_mode_list);
num_result_slots = numel(ilim_values) + compare_models.pn + num_m1pn_results;
if num_result_slots == 0
    error('m1_plane_source_benchmark:noModelsSelected', ...
        'At least one supported model must be enabled in compare_models.');
end

results = repmat(struct('name', '', 'rho', [], 'j', [], 'm2', [], ...
    'u0', [], 'u1', [], 'u2', [], ...
    'rho_l1', [], 'rho_linf', [], 'j_l1', [], 'j_linf', [], ...
    'm2_l1', [], 'm2_linf', [], 'diag', struct()), 1, num_result_slots);
result_idx = 0;

for k = 1:numel(ilim_values)
    ilim = ilim_values(k);
    u = simulate_m1(rho0, dt, dz, num_steps, ilim, sigma_a, sigma_s, source_strength, vacuum_opts);
    m2 = calc_psi2(u);

    result_idx = result_idx + 1;
    results(result_idx) = make_result(method_name(ilim), ...
        u(1, :)', u(2, :)', m2(:), ...
        u(1, :)', u(2, :)', nan(num_cells, 1), ...
        rho_ref, j_ref, m2_ref, dz, struct());
end

if compare_models.pn
    uPN = simulate_pn(rho0, dt, dz, num_steps, pn_order, pn_ilim, ...
        sigma_a, sigma_s, source_strength, vacuum_opts);
    obs_pn = pn_observables(uPN);
    result_idx = result_idx + 1;
    results(result_idx) = make_result(sprintf('PN-%s', method_name(pn_ilim)), ...
        obs_pn.rho, obs_pn.j, obs_pn.m2, ...
        uPN(1, :)', uPN(2, :)', pn_state_component(uPN, 3), ...
        rho_ref, j_ref, m2_ref, dz, struct());
end

if compare_models.m1pn
    m1pn_pn_ho_method_id = normalize_flux_method(m1pn_pn_ho_method);
    m1pn_pn_ho_method_name = method_name(m1pn_pn_ho_method_id);
    for mode_idx = 1:numel(m1pn_mode_list)
        m1pn_mode = m1pn_mode_list{mode_idx};
        m1pn_opts = vacuum_opts;
        m1pn_opts.pn_ho_method = m1pn_pn_ho_method_id;
        m1pn_opts.coupling_mode = m1pn_mode;
        m1pn_opts = merge_structs(m1pn_opts, m1pn_mode_opts);
        [u_m1pn, uPN_m1pn, m1pn_diag] = simulate_m1pn(rho0, dt, dz, num_steps, m1pn_order, ...
            sigma_a, sigma_s, source_strength, m1pn_opts);
        m2_m1pn = m1pn_second_moment(uPN_m1pn);
        result_idx = result_idx + 1;
        results(result_idx) = make_result(m1pn_result_name(m1pn_mode, m1pn_pn_ho_method_name), ...
            u_m1pn(1, :)', u_m1pn(2, :)', m2_m1pn(:), ...
            uPN_m1pn(1, :)', uPN_m1pn(2, :)', pn_state_component(uPN_m1pn, 3), ...
            rho_ref, j_ref, m2_ref, dz, summarize_m1pn_diag(m1pn_diag));
    end
end

results = results(1:result_idx);

%% Output
fprintf('Plane-source benchmark against master_basis S_N reference at T = %.3f\n', case_cfg.T);
for k = 1:numel(results)
    fprintf(['%s: rho L1 = %.6e, rho Linf = %.6e, ' ...
        'j L1 = %.6e, j Linf = %.6e, m2 L1 = %.6e, m2 Linf = %.6e\n'], ...
        results(k).name, ...
        results(k).rho_l1, results(k).rho_linf, ...
        results(k).j_l1, results(k).j_linf, ...
        results(k).m2_l1, results(k).m2_linf);
    if ~isempty(fieldnames(results(k).diag))
        fprintf(['  diag: fallback steps = %d, mean outer iters = %.3f, ' ...
            'max sync correction = %.6e\n'], ...
            results(k).diag.fallback_steps, ...
            results(k).diag.mean_outer_iterations, ...
            results(k).diag.max_sync_correction);
    end
end

%% Plots
plot_mask = z >= 0;

figure;
subplot(3, 1, 1)
hold on
plot(z(plot_mask), rho_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'S_N ref')
for k = 1:numel(results)
    plot(z(plot_mask), results(k).rho(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: \rho = <u>')
xlabel('z')

subplot(3, 1, 2)
hold on
plot(z(plot_mask), j_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'S_N ref')
for k = 1:numel(results)
    plot(z(plot_mask), results(k).j(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: j = <\mu u>')
xlabel('z')

subplot(3, 1, 3)
hold on
plot(z(plot_mask), m2_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'S_N ref')
for k = 1:numel(results)
    plot(z(plot_mask), results(k).m2(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: m_2 = <\mu^2 u>')
xlabel('z')

figure;
subplot(3, 1, 1)
hold on
plot(z(plot_mask), rho_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'u_0 ref')
for k = 1:numel(results)
    plot(z(plot_mask), results(k).u0(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: moment u_0')
xlabel('z')

subplot(3, 1, 2)
hold on
plot(z(plot_mask), j_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'u_1 ref')
for k = 1:numel(results)
    plot(z(plot_mask), results(k).u1(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: moment u_1')
xlabel('z')

subplot(3, 1, 3)
hold on
plot(z(plot_mask), u2_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'u_2 ref')
for k = 1:numel(results)
    if all(isnan(results(k).u2))
        continue;
    end
    plot(z(plot_mask), results(k).u2(plot_mask), 'DisplayName', results(k).name);
end
hold off
legend('Location', 'best')
title('Plane source: moment u_2 (PN state)')
xlabel('z')

function name = method_name(ilim)
    switch ilim
        case 1
            name = 'MCL';
        case -1
            name = 'LLF';
        otherwise
            name = 'LW';
    end
end

function name = m1pn_result_name(coupling_mode, pn_method_name)
    switch lower(strtrim(coupling_mode))
        case 'explicit_failsafe'
            coupling_label = 'M1PN-HOLO-explicit';
        case 'explicit_synced'
            coupling_label = 'M1PN-HOLO-sync';
        case 'implicit_holo'
            coupling_label = 'M1PN-HOLO-implicit';
        otherwise
            coupling_label = sprintf('M1PN-%s', coupling_mode);
    end
    name = sprintf('%s/PN-%s', coupling_label, pn_method_name);
end

function mode_list = normalize_mode_list(raw_modes)
    if ischar(raw_modes) || isstring(raw_modes)
        mode_list = cellstr(raw_modes);
        return;
    end

    if ~iscell(raw_modes) || isempty(raw_modes)
        error('m1_plane_source_benchmark:invalidM1PNModes', ...
            'm1pn_coupling_modes must be a non-empty string or cell array of strings.');
    end

    mode_list = raw_modes;
end

function method_id = normalize_flux_method(raw_method)
    if isstring(raw_method)
        raw_method = char(raw_method);
    end

    if ischar(raw_method)
        switch upper(strtrim(raw_method))
            case 'MCL'
                method_id = 1;
            case 'LLF'
                method_id = -1;
            case 'LW'
                method_id = -2;
            otherwise
                error('m1_plane_source_benchmark:invalidM1PNPNMethod', ...
                    'Unknown m1pn_pn_ho_method "%s".', raw_method);
        end
        return;
    end

    if ~isscalar(raw_method)
        error('m1_plane_source_benchmark:invalidM1PNPNMethod', ...
            'm1pn_pn_ho_method must be a scalar or one of the strings MCL/LLF/LW.');
    end

    if any(raw_method == [1, -1, -2])
        method_id = raw_method;
    else
        error('m1_plane_source_benchmark:invalidM1PNPNMethod', ...
            'Unsupported numeric m1pn_pn_ho_method %g.', raw_method);
    end
end

function merged = merge_structs(base_struct, override_struct)
    merged = base_struct;
    if isempty(override_struct)
        return;
    end

    override_fields = fieldnames(override_struct);
    for k = 1:numel(override_fields)
        field_name = override_fields{k};
        merged.(field_name) = override_struct.(field_name);
    end
end

function diag_summary = summarize_m1pn_diag(diag_struct)
    diag_summary = struct();
    diag_summary.fallback_steps = nnz(diag_struct.fallback_used);
    diag_summary.mean_outer_iterations = mean(diag_struct.outer_iterations);
    diag_summary.max_sync_correction = max(diag_struct.sync_correction);
end

function result = make_result(name, rho, j, m2, u0, u1, u2, rho_ref, j_ref, m2_ref, dz, diag_summary)
    result = struct();
    result.name = name;
    result.rho = rho(:);
    result.j = j(:);
    result.m2 = m2(:);
    result.u0 = u0(:);
    result.u1 = u1(:);
    result.u2 = u2(:);
    result.rho_l1 = dz * sum(abs(result.rho - rho_ref));
    result.rho_linf = max(abs(result.rho - rho_ref));
    result.j_l1 = dz * sum(abs(result.j - j_ref));
    result.j_linf = max(abs(result.j - j_ref));
    result.m2_l1 = dz * sum(abs(result.m2 - m2_ref));
    result.m2_linf = max(abs(result.m2 - m2_ref));
    result.diag = diag_summary;
end

function component = pn_state_component(uPN, idx)
    if size(uPN, 1) < idx
        component = nan(size(uPN, 2), 1);
        return;
    end

    component = uPN(idx, :).';
end
