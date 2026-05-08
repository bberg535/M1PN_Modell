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
    'pn', true, ...
    'm1pn', true);

pn_ilim = -2;
pn_order = 9;
m1pn_order = 13;
m1pn_pn_ho_method = 1;

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

num_result_slots = numel(ilim_values) + compare_models.pn + compare_models.m1pn;
if num_result_slots == 0
    error('m1_plane_source_benchmark:noModelsSelected', ...
        'At least one supported model must be enabled in compare_models.');
end

results = repmat(struct('name', '', 'rho', [], 'j', [], 'm2', [], ...
    'rho_l1', [], 'rho_linf', [], 'j_l1', [], 'j_linf', [], ...
    'm2_l1', [], 'm2_linf', []), 1, num_result_slots);
result_idx = 0;

for k = 1:numel(ilim_values)
    ilim = ilim_values(k);
    u = simulate_m1(rho0, dt, dz, num_steps, ilim, sigma_a, sigma_s, source_strength, vacuum_opts);
    m2 = calc_psi2(u);

    result_idx = result_idx + 1;
    results(result_idx) = make_result(method_name(ilim), ...
        u(1, :)', u(2, :)', m2(:), ...
        rho_ref, j_ref, m2_ref, dz);
end

if compare_models.pn
    uPN = simulate_pn(rho0, dt, dz, num_steps, pn_order, pn_ilim, ...
        sigma_a, sigma_s, source_strength, vacuum_opts);
    obs_pn = pn_observables(uPN);
    result_idx = result_idx + 1;
    results(result_idx) = make_result(sprintf('PN-%s', method_name(pn_ilim)), ...
        obs_pn.rho, obs_pn.j, obs_pn.m2, ...
        rho_ref, j_ref, m2_ref, dz);
end

if compare_models.m1pn
    m1pn_opts = vacuum_opts;
    m1pn_opts.pn_ho_method = m1pn_pn_ho_method;
    [u_m1pn, ~] = simulate_m1pn(rho0, dt, dz, num_steps, m1pn_order, ...
        sigma_a, sigma_s, source_strength, m1pn_opts);
    m2_m1pn = calc_psi2(u_m1pn);
    result_idx = result_idx + 1;
    results(result_idx) = make_result('M1PN-HOLO', ...
        u_m1pn(1, :)', u_m1pn(2, :)', m2_m1pn(:), ...
        rho_ref, j_ref, m2_ref, dz);
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

function result = make_result(name, rho, j, m2, rho_ref, j_ref, m2_ref, dz)
    result = struct();
    result.name = name;
    result.rho = rho(:);
    result.j = j(:);
    result.m2 = m2(:);
    result.rho_l1 = dz * sum(abs(result.rho - rho_ref));
    result.rho_linf = max(abs(result.rho - rho_ref));
    result.j_l1 = dz * sum(abs(result.j - j_ref));
    result.j_linf = max(abs(result.j - j_ref));
    result.m2_l1 = dz * sum(abs(result.m2 - m2_ref));
    result.m2_linf = max(abs(result.m2 - m2_ref));
end
