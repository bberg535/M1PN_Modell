% M1 plane-source benchmark against a high-resolution S_N reference
% based on Section 6.1.1 of Seminarquelle 2.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
paths = setup_project_paths();

%% Benchmark setup
compare_models = struct( ...
    'mcl', false, ...
    'llf', false, ...
    'lw', false, ...
    'pn', true, ...
    'm1pn', false);

pn_ilim = -2;
m1pn_ilim = -1;
pn_order = 9;
m1pn_order = 13;

dt = 0.001;
dz = 0.002;
T = 1.0;
domain = [-1.2, 1.2];
psi_vac = 0.5e-8;

sn_ref_factor = 1;
sn_ref_n_mu = 128;
sn_ref_cfl = 0.45;

domain_length = domain(2) - domain(1);
num_cells = round(domain_length / dz);
if abs(num_cells * dz - domain_length) > 1e-12 * max(1, domain_length)
    error('m1_plane_source_benchmark:meshMismatch', ...
        'Domain length %.16g is not an integer multiple of dz=%.16g.', ...
        domain_length, dz);
end
if mod(num_cells, 2) ~= 0
    error('m1_plane_source_benchmark:oddCellCount', ...
        'Plane-source benchmark requires an even number of cells.');
end

num_steps = round(T / dt);
if abs(num_steps * dt - T) > 1e-12 * max(1, T)
    error('m1_plane_source_benchmark:timeMismatch', ...
        'T=%.16g is not an integer multiple of dt=%.16g.', T, dt);
end

z = linspace(domain(1) + 0.5 * dz, domain(2) - 0.5 * dz, num_cells);
rho0 = (2 * psi_vac) * ones(num_cells, 1);
rho0(num_cells / 2) = rho0(num_cells / 2) + 1.0 / (2.0 * dz);
rho0(num_cells / 2 + 1) = rho0(num_cells / 2 + 1) + 1.0 / (2.0 * dz);

flux = @(u) ([u(2,:); u(3,:)]);
fluxj = @(u) 1;
sigma_a = zeros(1, num_cells);
sigma_s = ones(1, num_cells);
source_strength = 0.0;

%% Reference solution
sn_ref_cfg = struct();
sn_ref_cfg.domain = domain;
sn_ref_cfg.tf = T;
sn_ref_cfg.n_cells = sn_ref_factor * num_cells;
sn_ref_cfg.n_mu = sn_ref_n_mu;
sn_ref_cfg.cfl = sn_ref_cfl;
sn_ref_cfg.sigma_a = 0.0;
sn_ref_cfg.sigma_s = 1.0;
sn_ref_cfg.psi_vac_density = psi_vac;
sn_ref = sn_reference_plane_source(sn_ref_cfg);

rho_ref = coarsen_uniform_field(sn_ref.rho, sn_ref_factor);
j_ref = coarsen_uniform_field(sn_ref.j, sn_ref_factor);
m2_ref = coarsen_uniform_field(sn_ref.m2, sn_ref_factor);

%% Compare transport fluxes
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

include_m1pn_result = compare_models.m1pn && (m1pn_ilim == -1);
num_result_slots = numel(ilim_values) + compare_models.pn + include_m1pn_result;
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
    method = select_flux_method(ilim);

    u = zeros(2, num_cells);
    u(1, :) = rho0.';

    for step = 1:num_steps
        u = advance_m1(u, dt, dz, flux, fluxj, num_cells, ...
            sigma_a, sigma_s, method, source_strength);

        if any(~isfinite(u(:))) || any(~isreal(u(:)))
            error('m1_plane_source_benchmark:invalidState', ...
                'Method %s produced an invalid state at step %d.', ...
                method_name(ilim), step);
        end
    end

    psi2 = calc_psi2(u);
    result_idx = result_idx + 1;
    results(result_idx) = make_result(method_name(ilim), ...
        u(1, :)', u(2, :)', psi2(:), ...
        rho_ref, j_ref, m2_ref, dz);
end

if compare_models.pn
    uPN = simulate_periodic_pn(rho0, dt, dz, num_steps, pn_order, pn_ilim, ...
        sigma_a, sigma_s, source_strength);
    obs_pn = pn_observables(uPN(:, 1:num_cells));
    result_idx = result_idx + 1;
    results(result_idx) = make_result(sprintf('PN-%s', method_name(pn_ilim)), ...
        obs_pn.rho, obs_pn.j, obs_pn.m2, ...
        rho_ref, j_ref, m2_ref, dz);
end

if compare_models.m1pn
    if m1pn_ilim == -1
        [u_m1pn, uPN_m1pn] = simulate_periodic_m1pn(rho0, dt, dz, num_steps, ...
            m1pn_order, m1pn_ilim, sigma_a, sigma_s, source_strength);
        obs_m1pn = pn_closure_observables(uPN_m1pn(:, 1:num_cells), ...
            u_m1pn(:, 1:num_cells));
        result_idx = result_idx + 1;
        results(result_idx) = make_result(sprintf('M1PN-%s', method_name(m1pn_ilim)), ...
            u_m1pn(1, 1:num_cells)', u_m1pn(2, 1:num_cells)', obs_m1pn.m2, ...
            rho_ref, j_ref, m2_ref, dz);
    else
        fprintf(['M1PN-%s: unsupported ' ...
            '(coupled M1PN currently supports only LLF / ilim = -1)\n'], ...
            method_name(m1pn_ilim));
    end
end

results = results(1:result_idx);

%% Output
fprintf('Plane-source benchmark against S_N reference at T = %.3f\n', T);
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

function coarse = coarsen_uniform_field(fine, factor)
    fine = fine(:);
    if mod(numel(fine), factor) ~= 0
        error('m1_plane_source_benchmark:coarsenMismatch', ...
            'Field length %d is not divisible by factor %d.', numel(fine), factor);
    end
    coarse = mean(reshape(fine, factor, []), 1).';
end

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
