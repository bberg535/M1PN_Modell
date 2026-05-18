% Coupled HOLO M1PN model for the linear Boltzmann equation
% with periodic boundary conditions.
%
% Active coupling logic:
% - M1 is the resolved model for (rho, j),
% - PN is the auxiliary high-order transport model,
% - the M1 interface flux is computed as usual and then its second
%   component is replaced by the PN-based u_2 approximation,
% - PN is projected back onto the M1 low moments after each step.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

%% Inputs
cfg = struct();
cfg.cont_plot = false;
cfg.compute_sn_reference = true;
cfg.sn_ref_n_mu = 128;
cfg.sn_ref_cfl = 0.45;

pn_order = 3;
pn_ho_method = 1;
dt = 0.001;
dz = 0.002;
T = 1.0;
nec = round(2 / dz);
z = 0:dz:2;

sigma_a = dz * ones(1, nec);
sigma_s = dz * ones(1, nec);
source_strength = 1.0;
rho0 = [1 / dz; zeros(nec - 1, 1)];

opts = struct('boundary', 'periodic', 'psi_boundary', 0.0, 'pn_ho_method', pn_ho_method);
u = zeros(2, nec);
uPN = zeros(pn_order + 1, nec);
u(1, :) = rho0.';
uPN(1, :) = rho0.';

if cfg.cont_plot
    figure;
    plot_state(z, u, m1pn_second_moment(uPN));
end

num_steps = round(T / dt);
if abs(num_steps * dt - T) > 1e-12 * max(1, T)
    error('m1pn_model:timeGridMismatch', ...
        'T=%.16g is not an integer multiple of dt=%.16g.', T, dt);
end

for step = 1:num_steps
    t = step * dt;
    step_info = struct('mode', 'M1PN', 'step', step, 'time', t);
    [u, uPN, step_diag] = heuns_method(u, uPN, dt, dz, sigma_a, sigma_s, ...
        pn_order, source_strength, opts, step_info);
    [uPN, sync_diag] = pn_sync_low_moments(uPN, u(1:2, :));

    mass = sum(u(1, :)) * dz;
    pn_diag = pn_observables(uPN);
    [min_pn_reconstruction, idx] = min(pn_diag.min_reconstruction);

    fprintf(['t = %.3f, mass = %.10e, alpha_min = %.3e, alpha_max = %.3e, ' ...
        'limited interfaces = %d, max sync correction = %.3e\n'], ...
        t, mass, step_diag.min_alpha, step_diag.max_alpha, step_diag.num_limited, ...
        sync_diag.max_correction);
    fprintf('min rho_M1 = %.3e, min PN reconstruction = %.3e at %d\n', ...
        min(u(1, :)), min_pn_reconstruction, idx);

    if cfg.cont_plot
        plot_state(z, u, m1pn_second_moment(uPN));
        drawnow;
    end
end

%% Diagnostics and plotting
m2_model = m1pn_second_moment(uPN);
u_plot = [u, u(:, 1)];
m2_plot = [m2_model, m2_model(1)];
u0_plot = periodic_component(uPN, 1);
u1_plot = periodic_component(uPN, 2);
u2_plot = periodic_component(uPN, 3);

if cfg.compute_sn_reference
    sn_ref_cfg = struct();
    sn_ref_cfg.domain = [0.0, 2.0];
    sn_ref_cfg.tf = T;
    sn_ref_cfg.n_cells = nec;
    sn_ref_cfg.n_mu = cfg.sn_ref_n_mu;
    sn_ref_cfg.cfl = cfg.sn_ref_cfl;
    sn_ref_cfg.sigma_a = sigma_a;
    sn_ref_cfg.sigma_s = sigma_s;
    sn_ref_cfg.source_strength = source_strength;
    sn_ref_cfg.initial_rho = rho0(1:nec).';
    sn_ref = sn_reference_periodic_source(sn_ref_cfg);
    sn_ref_moments = sn_reference_pn_moments(sn_ref, 2);

    rho_ref = [sn_ref.rho; sn_ref.rho(1)];
    j_ref = [sn_ref.j; sn_ref.j(1)];
    m2_ref = [sn_ref.m2; sn_ref.m2(1)];
    u2_ref = [sn_ref_moments(3, :), sn_ref_moments(3, 1)];

    rho_err = abs(u(1, :).'- sn_ref.rho);
    j_err = abs(u(2, :).'- sn_ref.j);
    m2_err = abs(m2_model.' - sn_ref.m2);

    fprintf(['SN reference errors: rho L1 = %.6e, rho Linf = %.6e, ' ...
        'j L1 = %.6e, j Linf = %.6e, m2 L1 = %.6e, m2 Linf = %.6e\n'], ...
        dz * sum(rho_err), max(rho_err), ...
        dz * sum(j_err), max(j_err), ...
        dz * sum(m2_err), max(m2_err));
end

figure;
hold on
plot(z, u_plot(1, :), 'DisplayName', 'rho M1')
plot(z, u_plot(2, :), 'DisplayName', 'j M1')
plot(z, m2_plot, 'DisplayName', 'm_2 PN-projected')
if cfg.compute_sn_reference
    plot(z, rho_ref, '--', 'DisplayName', 'rho S_N')
    plot(z, j_ref, '--', 'DisplayName', 'j S_N')
    plot(z, m2_ref, '--', 'DisplayName', 'm_2 S_N')
end
hold off
legend('Location', 'best')

figure;
subplot(3, 1, 1)
hold on
plot(z, u0_plot, 'DisplayName', 'u_0 PN')
plot(z, u_plot(1, :), '--', 'DisplayName', 'u_0 M1')
if cfg.compute_sn_reference
    plot(z, rho_ref, ':', 'DisplayName', 'u_0 S_N')
end
hold off
legend('Location', 'best')
title('Moment u_0')

subplot(3, 1, 2)
hold on
plot(z, u1_plot, 'DisplayName', 'u_1 PN')
plot(z, u_plot(2, :), '--', 'DisplayName', 'u_1 M1')
if cfg.compute_sn_reference
    plot(z, j_ref, ':', 'DisplayName', 'u_1 S_N')
end
hold off
legend('Location', 'best')
title('Moment u_1')

subplot(3, 1, 3)
hold on
plot(z, u2_plot, 'DisplayName', 'u_2 PN')
plot(z, m2_plot, '--', 'DisplayName', 'm_2 PN-projected')
if cfg.compute_sn_reference
    plot(z, u2_ref, ':', 'DisplayName', 'u_2 S_N->PN')
end
hold off
legend('Location', 'best')
title('Moment u_2')

function plot_state(z, u, m2)
    state_plot = [u, u(:, 1)];
    m2_plot_local = [m2, m2(1)];
    plot(z, [state_plot; m2_plot_local])
end

function values = periodic_component(state, idx)
    if size(state, 1) < idx
        values = nan(1, size(state, 2) + 1);
        return;
    end

    values = [state(idx, :), state(idx, 1)];
end
