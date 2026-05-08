% PN vs M1PN negativity benchmark on the plane-source test.
%
% This uses the 1D slab setup from Section 6.1.1 ("plane source") of
% Seminarquelle 2. That section explicitly discusses strong oscillations
% for PN-type models on this non-smooth benchmark. The paper also reports
% negative PN particle densities in a later multi-D test; since this code
% base is 1D, the plane-source case is the closest directly reusable
% stress test here.
%
% Important: the current coupled M1PN implementation supports only LLF.
% To expose the positivity advantage of the coupling, this script runs
% pure PN with a low-diffusion LW flux and compares it against M1PN-LLF.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

%% Switches
include_sn_reference = true;
pn_ilim = -2;
m1pn_ilim = -1;
pn_order = 9;
m1pn_order = pn_order;
assert_negative_pn = true;
assert_positive_m1pn = true;

%% Plane-source setup
dt = 1.0e-3;
dz = 4.0e-3;
T = 1.0;
domain = [-1.2, 1.2];
psi_vac = 0.5e-8;
sigma_a_value = 0.0;
sigma_s_value = 1.0;
source_strength = 0.0;

sn_ref_factor = 2;
sn_ref_n_mu = 128;
sn_ref_cfl = 0.45;
zoom_half_width = 0.2;

if m1pn_ilim ~= -1
    error('pn_m1pn_negativity_benchmark:unsupportedCoupledMethod', ...
        'M1PN currently supports only ilim = -1 (LLF).');
end

[z, rho0, num_cells, num_steps, sigma_a, sigma_s] = plane_source_case( ...
    domain, dz, T, dt, psi_vac, sigma_a_value, sigma_s_value);

%% Pure PN
uPN = simulate_periodic_pn(rho0, dt, dz, num_steps, pn_order, pn_ilim, ...
    sigma_a, sigma_s, source_strength);
obs_pn = pn_observables(uPN(:, 1:num_cells));

%% Coupled M1PN
[u_m1pn, uPN_m1pn] = simulate_periodic_m1pn(rho0, dt, dz, num_steps, ...
    m1pn_order, m1pn_ilim, sigma_a, sigma_s, source_strength);
obs_m1pn = pn_closure_observables(uPN_m1pn(:, 1:num_cells), ...
    u_m1pn(:, 1:num_cells));
rho_m1pn = u_m1pn(1, 1:num_cells).';

%% Optional S_N reference
if include_sn_reference
    sn_ref_cfg = struct();
    sn_ref_cfg.domain = domain;
    sn_ref_cfg.tf = T;
    sn_ref_cfg.n_cells = sn_ref_factor * num_cells;
    sn_ref_cfg.n_mu = sn_ref_n_mu;
    sn_ref_cfg.cfl = sn_ref_cfl;
    sn_ref_cfg.sigma_a = sigma_a_value;
    sn_ref_cfg.sigma_s = sigma_s_value;
    sn_ref_cfg.psi_vac_density = psi_vac;

    sn_ref = sn_reference_plane_source(sn_ref_cfg);
    rho_ref = coarsen_uniform_field(sn_ref.rho, sn_ref_factor);
end

%% Diagnostics
[pn_min_rho, pn_min_idx] = min(obs_pn.rho);
[m1pn_min_rho, m1pn_min_idx] = min(rho_m1pn);

pn_negative_cells = nnz(obs_pn.rho < 0);
m1pn_negative_cells = nnz(rho_m1pn < 0);

fprintf('PN vs M1PN negativity benchmark (plane source, Section 6.1.1)\n');
fprintf('PN-%s (N=%d):   min rho = %.6e at z = %.6f, negative cells = %d\n', ...
    method_name(pn_ilim), pn_order, pn_min_rho, z(pn_min_idx), pn_negative_cells);
fprintf('M1PN-%s (N=%d): min rho = %.6e at z = %.6f, negative cells = %d\n', ...
    method_name(m1pn_ilim), m1pn_order, m1pn_min_rho, z(m1pn_min_idx), ...
    m1pn_negative_cells);

if include_sn_reference
    fprintf('PN-%s rho L1 vs S_N   = %.6e\n', ...
        method_name(pn_ilim), dz * sum(abs(obs_pn.rho - rho_ref)));
    fprintf('M1PN-%s rho L1 vs S_N = %.6e\n', ...
        method_name(m1pn_ilim), dz * sum(abs(rho_m1pn - rho_ref)));
end

if assert_negative_pn && pn_min_rho >= 0
    error('pn_m1pn_negativity_benchmark:missingPNNegativity', ...
        'Expected negative PN density, but min(rho_PN) = %.16e.', pn_min_rho);
end

if assert_positive_m1pn && m1pn_min_rho < 0
    error('pn_m1pn_negativity_benchmark:negativeM1PNDensity', ...
        'Expected non-negative M1PN density, but min(rho_M1PN) = %.16e.', ...
        m1pn_min_rho);
end

%% Plots
plot_mask = z >= 0;
zoom_center = z(pn_min_idx);
zoom_mask = plot_mask & (abs(z - zoom_center) <= zoom_half_width);
if nnz(zoom_mask) < 5
    zoom_mask = plot_mask;
end

figure;

subplot(2, 1, 1)
hold on
if include_sn_reference
    plot(z(plot_mask), rho_ref(plot_mask), 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'S_N ref');
end
plot(z(plot_mask), obs_pn.rho(plot_mask), 'Color', [0.85, 0.33, 0.10], ...
    'LineWidth', 1.1, ...
    'DisplayName', sprintf('PN-%s', method_name(pn_ilim)));
plot(z(plot_mask), rho_m1pn(plot_mask), 'Color', [0.00, 0.45, 0.74], ...
    'LineWidth', 1.1, ...
    'DisplayName', sprintf('M1PN-%s', method_name(m1pn_ilim)));
yline(0.0, 'k--', 'DisplayName', 'rho = 0');
hold off
legend('Location', 'best')
title('Plane source: rho = <u> (right half-domain)')
xlabel('z')
ylabel('rho')

subplot(2, 1, 2)
hold on
if include_sn_reference
    plot(z(zoom_mask), rho_ref(zoom_mask), 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'S_N ref');
end
plot(z(zoom_mask), obs_pn.rho(zoom_mask), 'Color', [0.85, 0.33, 0.10], ...
    'LineWidth', 1.1, 'DisplayName', sprintf('PN-%s', method_name(pn_ilim)));
plot(z(zoom_mask), rho_m1pn(zoom_mask), 'Color', [0.00, 0.45, 0.74], ...
    'LineWidth', 1.1, 'DisplayName', sprintf('M1PN-%s', method_name(m1pn_ilim)));
yline(0.0, 'k--', 'DisplayName', 'rho = 0');
hold off
legend('Location', 'best')
title(sprintf('Zoom around PN minimum at z = %.3f', zoom_center))
xlabel('z')
ylabel('rho')

function [z, rho0, num_cells, num_steps, sigma_a, sigma_s] = plane_source_case( ...
    domain, dz, T, dt, psi_vac, sigma_a_value, sigma_s_value)
    domain_length = domain(2) - domain(1);
    num_cells = round(domain_length / dz);
    if abs(num_cells * dz - domain_length) > 1e-12 * max(1, domain_length)
        error('pn_m1pn_negativity_benchmark:meshMismatch', ...
            'Domain length %.16g is not an integer multiple of dz=%.16g.', ...
            domain_length, dz);
    end
    if mod(num_cells, 2) ~= 0
        error('pn_m1pn_negativity_benchmark:oddCellCount', ...
            'Plane-source benchmark requires an even number of cells.');
    end

    num_steps = round(T / dt);
    if abs(num_steps * dt - T) > 1e-12 * max(1, T)
        error('pn_m1pn_negativity_benchmark:timeMismatch', ...
            'T=%.16g is not an integer multiple of dt=%.16g.', T, dt);
    end

    boundary_distance = min(abs(domain));
    if T >= boundary_distance
        error('pn_m1pn_negativity_benchmark:boundaryReach', ...
            ['This benchmark reuses the periodic PN/M1PN stepping under ' ...
            'the assumption T < distance to the boundary.']);
    end

    z = linspace(domain(1) + 0.5 * dz, domain(2) - 0.5 * dz, num_cells).';
    rho0 = (2 * psi_vac) * ones(num_cells, 1);
    rho0(num_cells / 2) = rho0(num_cells / 2) + 1.0 / (2.0 * dz);
    rho0(num_cells / 2 + 1) = rho0(num_cells / 2 + 1) + 1.0 / (2.0 * dz);

    sigma_a = sigma_a_value * ones(1, num_cells);
    sigma_s = sigma_s_value * ones(1, num_cells);
end

function coarse = coarsen_uniform_field(fine, factor)
    fine = fine(:);
    if mod(numel(fine), factor) ~= 0
        error('pn_m1pn_negativity_benchmark:coarsenMismatch', ...
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
