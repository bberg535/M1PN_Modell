function run_planesource_fig1()
% Reproduction attempt for Figure 1 (plane-source test, slab geometry, tf=1)
% from "Seminarquelle 2".
%
% Produces rho(x) for MN and PN with N in {3,7,51} and a high-res SN reference.

%clc; close all;

%% Parameters
tf      = 1.0;
L       = 1.2;                % domain [-L, L] ;
Nx      = 400;
cfl     = 0.9;

psi_vac = 0.5e-8;             % vacuum density
rho_vac = 2 * psi_vac;        % rho = <psi> on mu in [-1,1]

sigma_s = 1.0;                % scattering
sigma_a = 0.0;                % absorption
Q0      = 0.0;                % isotropic source

dx      = 2 * L / Nx;
xC      = linspace(-L+dx/2, L-dx/2, Nx);

dt      = cfl*dx;             % characteristic speeds bounded by 1 in slab
Nt      = ceil(tf/dt);
dt      = tf/Nt;

%% --- Reference: high-res discrete ordinates SN (for plotting "Ref")
ref.Nmu     = 240;            % increase if you want sharper reference
ref.dt      = dt;             % keep same dt for comparability
ref.Nx      = Nx;             % you can also set larger for reference
ref.dx      = dx;
ref.L       = L;
ref.psi_vac = psi_vac;
ref.sigma_s = sigma_s;
ref.sigma_a = sigma_a;
ref.Q0      = Q0;
rho_ref = solve_SN_reference(ref, tf);

%% --- Models/orders shown in Figure 1
Ns = [3];

rho_PN = cell(numel(Ns),1);
rho_MN = cell(numel(Ns),1);

for k = 1:numel(Ns)
    N = Ns(k);

    % Quadrature for moments/flux integrals
    qOrder = max(80, 2*N + 40);     % paper uses 2N+40 (common choice there)
    quad = make_quadrature(qOrder);

    % Solve PN
    opts = struct();
    opts.model   = 'PN';
    opts.N       = N;
    opts.L       = L;
    opts.Nx      = Nx;
    opts.dx      = dx;
    opts.dt      = dt;
    opts.Nt      = Nt;
    opts.psi_vac = psi_vac;
    opts.rho_vac = rho_vac;
    opts.sigma_s = sigma_s;
    opts.sigma_a = sigma_a;
    opts.Q0      = Q0;
    opts.quad    = quad;
    rho_PN{k} = solve_moment_model(opts);

    % Solve MN
    opts.model   = 'MN';
    opts.newton_maxit = 50;
    opts.newton_tol   = 1e-10;
    opts.newton_damp  = 0.5;
    rho_MN{k} = solve_moment_model(opts);
end

%% --- Plot like Figure 1: show only x in [0,1]
idx = (xC >= 0);

figure('Color','w','Position',[100 100 1100 420]);

subplot(1,2,1); hold on; box on;
plot(xC(idx), rho_ref(idx), 'k-', 'LineWidth', 1.2);
for k=1:numel(Ns)
    plot(xC(idx), rho_MN{k}(idx), 'LineWidth', 1.2);
end
title('M_N (minimum entropy)'); xlabel('x'); ylabel('\rho(x,t_f)');
legend_str = [{'Ref (S_N)'} arrayfun(@(N) sprintf('M_%d',N), Ns, 'UniformOutput',false)];
legend(legend_str{:}, 'Location','northeast');
xlim([0 1]);

subplot(1,2,2); hold on; box on;
plot(xC(idx), rho_ref(idx), 'k-', 'LineWidth', 1.2);
for k=1:numel(Ns)
    plot(xC(idx), rho_PN{k}(idx), 'LineWidth', 1.2);
end
title('P_N (spherical harmonics / Legendre truncation)'); xlabel('x'); ylabel('\rho(x,t_f)');
legend_str = [{'Ref (S_N)'} arrayfun(@(N) sprintf('P_%d',N), Ns, 'UniformOutput',false)];
legend(legend_str{:}, 'Location','northeast');
xlim([0 1]);

end

%% ===================== Core solvers =====================







%% ===================== PN closure (linear) =====================


%% ===================== MN closure (entropy, exp ansatz) =====================






%% ===================== SN reference solver =====================





%% ===================== Utilities =====================




