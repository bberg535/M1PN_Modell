% Apply the active HOLO/M1PN model to the Figure 3 test case from:
% C. Hauck and R. McClarren, "Positive PN Closures",
% SIAM J. Sci. Comput. 32(5), 2010, DOI 10.1137/090764918.
%
% The paper/legacy SIU-PN run is kept as a like-for-like reference on the
% same periodic plane-source setup. The actual comparison target here is
% the current M1PN model in this repository.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

cfg = struct();
cfg.N = 1;
cfg.dt = 1.0e-3;
cfg.dz = 2.0e-3;
cfg.tf = 1.0;
cfg.domain_length = 2.0;
cfg.sigma = 1.0;
cfg.mu_eval = linspace(-1.0, 1.0, 1000);
pn_method = 1;
m1_method = -1;
pn_ho_method = 1;

pn_label = 'PN-MCL';
m1_label = 'M1-LLF';

num_cells = round(cfg.domain_length / cfg.dz);
num_steps = round(cfg.tf / cfg.dt);
rho0 = zeros(1, num_cells);
rho0(1) = 1.0 / cfg.dz;

% Figure 3 uses pure scattering with periodic boundaries and no source.
sigma_a = 0.0;
sigma_s = cfg.sigma;
source_strength = 0.0;
opts = struct('boundary', 'periodic', 'psi_boundary', 0.0, ...
    'pn_ho_method', pn_ho_method);

pn_ref = simulate_pn_siu(cfg);
u_pn_mcl = simulate_pn(rho0, cfg.dt, cfg.dz, num_steps, cfg.N, ...
    pn_method, sigma_a, sigma_s, 0.0, opts);
u_m1 = simulate_m1(rho0, cfg.dt, cfg.dz, num_steps, m1_method, ...
    sigma_a, sigma_s, source_strength, opts);
[u_m1pn, uPN_aux, holo_diag] = simulate_m1pn(rho0, cfg.dt, cfg.dz, ...
    num_steps, cfg.N, sigma_a, sigma_s, source_strength, opts);

pn_ref_obs = pn_observables(pn_ref.u);
pn_mcl_obs = pn_observables(u_pn_mcl);

G = diag(2 ./ (2 .* (0:cfg.N) + 1));
p_eval = legpoly_eval(cfg.mu_eval, cfg.N);
pn_mcl_coeff = G \ u_pn_mcl;
pn_mcl_reconstruction = pn_mcl_coeff.' * p_eval;
aux_coeff = G \ uPN_aux;
aux_reconstruction = aux_coeff.' * p_eval;
[m1_reconstruction, m1_min_reconstruction] = m1_entropy_reconstruction(u_m1, cfg.mu_eval);
[m1pn_reconstruction, m1pn_min_reconstruction] = m1_entropy_reconstruction(u_m1pn, cfg.mu_eval);

pn_mcl_min_reconstruction = min(pn_mcl_reconstruction, [], 2);
pn_aux_min_reconstruction = min(aux_reconstruction, [], 2);
[~, pn_mcl_worst_idx] = min(pn_mcl_min_reconstruction);
[~, m1_worst_idx] = min(m1_min_reconstruction);
[~, m1pn_worst_idx] = min(m1pn_min_reconstruction);

paper_peak_profile_ref = pn_ref.worst_profile;
paper_peak_profile_pn_mcl = pn_mcl_reconstruction(pn_mcl_worst_idx, :).';
paper_peak_profile_m1 = m1_reconstruction(m1_worst_idx, :).';
paper_peak_profile_m1pn = m1pn_reconstruction(m1pn_worst_idx, :).';

rho_m1 = u_m1(1, :).';
j_m1 = u_m1(2, :).';
rho_pn_mcl = u_pn_mcl(1, :).';
j_pn_mcl = u_pn_mcl(2, :).';
rho_m1pn = u_m1pn(1, :).';
j_m1pn = u_m1pn(2, :).';
rho_l1_error_pn_mcl = cfg.dz * sum(abs(rho_pn_mcl - pn_ref.rho));
j_l1_error_pn_mcl = cfg.dz * sum(abs(j_pn_mcl - pn_ref_obs.j));
rho_l1_error_m1 = cfg.dz * sum(abs(rho_m1 - pn_ref.rho));
j_l1_error_m1 = cfg.dz * sum(abs(j_m1 - pn_ref_obs.j));
rho_l1_error = cfg.dz * sum(abs(rho_m1pn - pn_ref.rho));
j_l1_error = cfg.dz * sum(abs(j_m1pn - pn_ref_obs.j));
num_limited_steps = nnz(holo_diag.num_limited > 0);

fprintf('Positive PN Closures Figure 3 test case with active M1PN model\n');
fprintf('N = %d, dt = %.6g, dz = %.6g, tf = %.6g, sigma = %.6g\n', ...
    pn_ref.N, pn_ref.dt, pn_ref.dz, pn_ref.tf, pn_ref.sigma);
fprintf('Reference PN-SIU: min F* = %.6e at z = %.6f\n', ...
    min(pn_ref.min_reconstruction), pn_ref.worst_z);
fprintf('%s: rho L1 = %.6e, j L1 = %.6e, min F* = %.6e\n', ...
    pn_label, rho_l1_error_pn_mcl, j_l1_error_pn_mcl, ...
    min(pn_mcl_min_reconstruction));
fprintf('%s: rho L1 = %.6e, j L1 = %.6e, min psi_ME = %.6e\n', ...
    m1_label, rho_l1_error_m1, j_l1_error_m1, ...
    min(m1_min_reconstruction));
fprintf('M1PN: rho L1 = %.6e, j L1 = %.6e, min psi_ME = %.6e\n', ...
    rho_l1_error, j_l1_error, min(m1pn_min_reconstruction));
fprintf('M1PN auxiliary PN: min F* = %.6e\n', min(pn_aux_min_reconstruction));
fprintf(['M1PN limiter: min alpha = %.6e, max alpha = %.6e, ' ...
    'limited steps = %d / %d\n'], ...
    min(holo_diag.min_alpha), max(holo_diag.max_alpha), ...
    num_limited_steps, num_steps);
fprintf(['Panel (d) sample points: PN-SIU z = %.6f, PN-MCL z = %.6f, ' ...
    'M1 z = %.6f, M1PN z = %.6f\n'], ...
    pn_ref.worst_z, pn_ref.z(pn_mcl_worst_idx), ...
    pn_ref.z(m1_worst_idx), pn_ref.z(m1pn_worst_idx));

figure;

subplot(2, 2, 1)
hold on
plot(pn_ref.z, pn_ref.rho, 'k', 'LineWidth', 1.1, ...
    'DisplayName', 'PN-SIU reference')
plot(pn_ref.z, rho_pn_mcl, 'b-.', 'LineWidth', 1.1, ...
    'DisplayName', pn_label)
plot(pn_ref.z, rho_m1, 'g:', 'LineWidth', 1.2, ...
    'DisplayName', m1_label)
plot(pn_ref.z, rho_m1pn, 'r--', 'LineWidth', 1.1, ...
    'DisplayName', 'M1PN')
hold off
title('(a) Scalar flux, \rho')
xlabel('x')
legend('Location', 'best')

subplot(2, 2, 2)
hold on
plot(pn_ref.z, pn_ref_obs.j, 'k', 'LineWidth', 1.1, ...
    'DisplayName', 'PN-SIU reference')
plot(pn_ref.z, j_pn_mcl, 'b-.', 'LineWidth', 1.1, ...
    'DisplayName', pn_label)
plot(pn_ref.z, j_m1, 'g:', 'LineWidth', 1.2, ...
    'DisplayName', m1_label)
plot(pn_ref.z, j_m1pn, 'r--', 'LineWidth', 1.1, ...
    'DisplayName', 'M1PN')
hold off
title('(b) Flux, j = <\mu u>')
xlabel('x')
legend('Location', 'best')

subplot(2, 2, 3)
hold on
plot(pn_ref.z, pn_ref.min_reconstruction, 'k', 'LineWidth', 1.1, ...
    'DisplayName', 'PN-SIU reference')
plot(pn_ref.z, pn_mcl_min_reconstruction, 'b-.', 'LineWidth', 1.1, ...
    'DisplayName', pn_label)
plot(pn_ref.z, m1_min_reconstruction, 'g:', 'LineWidth', 1.2, ...
    'DisplayName', m1_label)
plot(pn_ref.z, m1pn_min_reconstruction, 'r--', 'LineWidth', 1.1, ...
    'DisplayName', 'M1PN')
plot(pn_ref.z, zeros(size(pn_ref.z)), 'k:', 'HandleVisibility', 'off')
hold off
title('(c) min_{\mu \in [-1,1]} \psi(\cdot,\mu)')
xlabel('x')
legend('Location', 'best')

subplot(2, 2, 4)
hold on
plot(pn_ref.mu_eval, paper_peak_profile_ref, 'k', 'LineWidth', 1.1, ...
    'DisplayName', 'PN-SIU reference')
plot(cfg.mu_eval, paper_peak_profile_pn_mcl, 'b-.', 'LineWidth', 1.1, ...
    'DisplayName', pn_label)
plot(cfg.mu_eval, paper_peak_profile_m1, 'g:', 'LineWidth', 1.2, ...
    'DisplayName', m1_label)
plot(cfg.mu_eval, paper_peak_profile_m1pn, 'r--', 'LineWidth', 1.1, ...
    'DisplayName', 'M1PN')
scatter(pn_ref.lambda, pn_ref.quadrature_values, 30, 'k', 'filled', ...
    'DisplayName', 'PN quadrature values')
scatter(pn_ref.lambda, pn_mcl_obs.angular_values(:, pn_mcl_worst_idx), 30, 'b', ...
    'filled', 'DisplayName', 'PN-MCL values')
plot(pn_ref.mu_eval, zeros(size(pn_ref.mu_eval)), 'k:', 'HandleVisibility', 'off')
hold off
title(sprintf('(d) Worst reconstruction profiles'))
xlabel('\mu')
legend('Location', 'best')
