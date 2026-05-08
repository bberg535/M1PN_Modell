% M1PN-Modell fuer die lineare Boltzmann-Gleichung
% mit periodischen Randwerten.

clear;
script_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(script_dir));
setup_project_paths();

%% Eingabewerte
% Verfahren bestimmen: -1 [LLF] ist der einzige unterstuetzte gekoppelte Pfad.
ilim = -1;
cont_plot = false;
compute_sn_reference = true;
sn_ref_n_mu = 128;
sn_ref_cfl = 0.45;

% Groessen fuer die Diskretisierung
NPN = 3;
dt = 0.001;
dz = 0.002;    % CFL: dt < dz
nec = 2 / dz;
nne = nec + 1;
T = 1;

z = 0:dz:2;

% Flussfunktion
flux = @(u) ([u(2,:); u(3,:)]);

% Absorptions- und Streuungsterm
sigma_a = (dz .* ones(nne, 1))';
sigma_s = (dz .* ones(nne, 1))';
source_strength = 1.0;

% Anfangswerte
rho0 = [1 / dz; zeros(nne - 1, 1)];
u = zeros(2, nne);
uPN = zeros(NPN + 1, nne);
u(1, :) = rho0.';
uPN(1, :) = rho0.';
u(:, nne) = u(:, 1);
uPN(:, nne) = uPN(:, 1);

if cont_plot
    figure;
    plot(z, [u; calc_psi2(u)])
end

%% Zeitschleife
num_steps = round(T / dt);
if abs(num_steps * dt - T) > 1e-12 * max(1, T)
    error('m1pn_model:timeGridMismatch', ...
        'T=%.16g is not an integer multiple of dt=%.16g.', T, dt);
end

for step = 1:num_steps
    t = step * dt;
    step_info = struct('mode', 'M1PN', 'step', step, 'time', t);

    % Zeitschritt mit SSP Heuns Method
    [u, uPN] = heuns_method(u, uPN, dt, dz, flux, nec, sigma_a, sigma_s, ...
        ilim, NPN, source_strength, step_info);

    mass = sum(u(1, 1:nec)) * dz;
    pn_diag = pn_observables(uPN(:, 1:nec));
    [min_pn_reconstruction, idx] = min(pn_diag.min_reconstruction);

    fprintf('t = %.3f, mass = %.10e\n', t, mass);
    fprintf('min rho_M1 = %.3e, min PN reconstruction = %.3e at %d\n ', ...
        min(u(1, 1:nec)), min_pn_reconstruction, idx);

    if cont_plot
        plot(z, [u; calc_psi2(u)])
        drawnow;
    end
end

%% Plotten
pn_obs = pn_observables(uPN(:, 1:nec));
closure_obs = pn_closure_observables(uPN(:, 1:nec), u(:, 1:nec));
m2_model = [closure_obs.m2.', closure_obs.m2(1)];

if compute_sn_reference
    sn_ref_cfg = struct();
    sn_ref_cfg.domain = [0.0, 2.0];
    sn_ref_cfg.tf = T;
    sn_ref_cfg.n_cells = nec;
    sn_ref_cfg.n_mu = sn_ref_n_mu;
    sn_ref_cfg.cfl = sn_ref_cfl;
    sn_ref_cfg.sigma_a = sigma_a(1:nec);
    sn_ref_cfg.sigma_s = sigma_s(1:nec);
    sn_ref_cfg.source_strength = source_strength;
    sn_ref_cfg.initial_rho = rho0(1:nec).';
    sn_ref = sn_reference_periodic_source(sn_ref_cfg);

    rho_ref = [sn_ref.rho; sn_ref.rho(1)];
    j_ref = [sn_ref.j; sn_ref.j(1)];
    m2_ref = [sn_ref.m2; sn_ref.m2(1)];

    rho_err = abs(u(1, 1:nec).' - sn_ref.rho);
    j_err = abs(u(2, 1:nec).' - sn_ref.j);
    m2_err = abs(closure_obs.m2 - sn_ref.m2);

    fprintf(['SN reference errors: rho L1 = %.6e, rho Linf = %.6e, ' ...
        'j L1 = %.6e, j Linf = %.6e, m2 L1 = %.6e, m2 Linf = %.6e\n'], ...
        dz * sum(rho_err), max(rho_err), ...
        dz * sum(j_err), max(j_err), ...
        dz * sum(m2_err), max(m2_err));
end

figure;
hold on
plot(z, u(1, :), 'DisplayName', 'rho M1')
plot(z, u(2, :), 'DisplayName', 'j M1')
plot(z, m2_model, 'DisplayName', 'm_2 PN closure')
if compute_sn_reference
    plot(z, rho_ref, '--', 'DisplayName', 'rho S_N')
    plot(z, j_ref, '--', 'DisplayName', 'j S_N')
    plot(z, m2_ref, '--', 'DisplayName', 'm_2 S_N')
end
hold off
legend('Location', 'best')

figure;
plot(z, uPN)
