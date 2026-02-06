%% run_planesource_fig1_MCL.m
% Plane-source test in slab geometry with PN and MN.
% Transport: MCL based on LLF (Rusanov) low-order and LW high-order flux.
%
% ilim =  1 : MCL (LLF + limited antidiffusion towards LW)
% ilim =  0 : pure LW
% ilim = -1 : pure LLF
%
% NOTE: For MN the MCL/LW transport can violate realizability; guard included.

clear; close all; clc;

%% ---------------- User parameters ----------------
tf      = 1.0;
L       = 1.2;            % domain [-L, L] (paper uses [-1.2,1.2] for tf=1)
Nx      = 400;           % even (Dirac split onto 2 cells)
cfl     = 0.5;            % dt/dx
ilim    = 1;              % 1=MCL, 0=LW, -1=LLF

psi_vac = 0.5e-8;
rho_vac = 2*psi_vac;

sigma_s = 1.0;
sigma_a = 0.0;
Q0      = 0.0;

Ns = [3];            % moments orders shown

% Quadrature for MN closure (only)
qPad = 40;                % use qOrder = 2N + qPad

%% ---------------- Grid/time ----------------
dx = 2*L/Nx;
xC = linspace(-L+dx/2, L-dx/2, Nx);

dt = cfl*dx;              % since max wave speed ~1
Nt = ceil(tf/dt);
dt = tf/Nt;
cfl = dt/dx;              % adjust to hit tf exactly

%% ---------------- Reference SN (optional) ----------------
doSNref = true;
rho_ref = [];
if doSNref
    ref.Nmu     = 240;
    ref.Nx      = Nx;
    ref.dx      = dx;
    ref.dt      = dt;
    ref.tf      = tf;
    ref.L       = L;
    ref.psi_vac = psi_vac;
    ref.sigma_s = sigma_s;
    ref.sigma_a = sigma_a;
    ref.Q0      = Q0;
    ref.ilim    = ilim;   % use same MCL for SN angles
    rho_ref = solve_SN_reference_MCL(ref);
end

%% ---------------- Solve PN/MN ----------------
rho_PN = cell(numel(Ns),1);
rho_MN = cell(numel(Ns),1);

for k=1:numel(Ns)
    N = Ns(k);

    % PN
    opts = struct();
    opts.model   = 'PN';
    opts.N       = N;
    opts.Nx      = Nx;
    opts.dx      = dx;
    opts.dt      = dt;
    opts.Nt      = Nt;
    opts.cfl     = cfl;
    opts.ilim    = ilim;
    opts.rho_vac = rho_vac;
    opts.psi_vac = psi_vac;
    opts.sigma_s = sigma_s;
    opts.sigma_a = sigma_a;
    opts.Q0      = Q0;

    rho_PN{k} = solve_moment_model_MCL(opts);

    % MN
    opts.model = 'MN';
    opts.qOrder = 2*N + qPad;
    opts.newton_maxit = 50;
    opts.newton_tol   = 1e-10;
    opts.newton_damp  = 0.5;

    rho_MN{k} = solve_moment_model_MCL(opts);
end

%% ---------------- Plot (x in [0,1]) ----------------
idx = (xC >= 0) & (xC <= 1);

figure('Color','w','Position',[100 100 1100 420]);

subplot(1,2,1); hold on; box on;
if doSNref
    plot(xC(idx), rho_ref(idx), 'k-', 'LineWidth', 1.2);
end
for k=1:numel(Ns)
    plot(xC(idx), rho_MN{k}(idx), 'LineWidth', 1.2);
end
title('M_N with MCL(LLF,LW) transport'); xlabel('x'); ylabel('\rho(x,t_f)');
leg = {};
if doSNref, leg{end+1}='Ref (S_N)'; end
for k=1:numel(Ns), leg{end+1}=sprintf('M_%d',Ns(k)); end
legend(leg{:}, 'Location','northeast'); xlim([0 1]);

subplot(1,2,2); hold on; box on;
if doSNref
    plot(xC(idx), rho_ref(idx), 'k-', 'LineWidth', 1.2);
end
for k=1:numel(Ns)
    plot(xC(idx), rho_PN{k}(idx), 'LineWidth', 1.2);
end
title('P_N with MCL(LLF,LW) transport'); xlabel('x'); ylabel('\rho(x,t_f)');
leg = {};
if doSNref, leg{end+1}='Ref (S_N)'; end
for k=1:numel(Ns), leg{end+1}=sprintf('P_%d',Ns(k)); end
legend(leg{:}, 'Location','northeast'); xlim([0 1]);

%% =======================================================================
%% ========================= Local functions =============================
%% =======================================================================

function rho = solve_moment_model_MCL(opts)
N   = opts.N;
nM  = N+1;
Nx  = opts.Nx;

% Initial moments: isotropic vacuum + Dirac split onto two middle cells
u = zeros(nM, Nx);
u(1,:) = opts.rho_vac;

midL = Nx/2; midR = Nx/2 + 1;
u(1,[midL midR]) = opts.rho_vac + 1/opts.dx;  % since rho = 2*(1/(2dx)) = 1/dx

% MN caches / quadrature
alpha_cache = [];
quad = [];
Pfull = [];
wmu_full = [];
if strcmpi(opts.model,'MN')
    quad = make_quadrature(opts.qOrder);
    Pfull = legendreP_matrix(N, quad.mu_full);
    wmu_full = quad.w_full(:).*quad.mu_full(:);
    alpha_cache = zeros(nM, Nx);
    alpha_cache(1,:) = log(max(u(1,:)/2, 1e-300)); % isotropic guess
end

for it=1:opts.Nt
    % half source step (exact scatter/absorption/isotropic Q)
    u = source_halfstep(u, opts.dt/2, opts.sigma_s, opts.sigma_a, opts.Q0);

    % transport step: MCL with LLF/LW flux
    if strcmpi(opts.model,'PN')
        u = transport_MCL_step(u, opts, @(uu) flux_PN(uu,N), [], []);
    else
        % update multipliers once per time step (cellwise Newton)
        [alpha_cache, ok] = multipliers_MN(u, alpha_cache, opts, Pfull, quad);
        if ~ok
            % guard: if Newton fails somewhere, we keep going (those cells get isotropized inside multipliers_MN)
        end
        f = flux_from_alpha(alpha_cache, Pfull, wmu_full);
        u = transport_MCL_step(u, opts, [], f, alpha_cache);
    end

    % half source step
    u = source_halfstep(u, opts.dt/2, opts.sigma_s, opts.sigma_a, opts.Q0);

    % realizability/positivity guard on density (helps MN stability)
    u(1,:) = max(u(1,:), opts.rho_vac);
end

rho = u(1,:);
end

function u_new = transport_MCL_step(u, opts, fluxfun_cell, f_cell, alpha_cache)
% One explicit transport step using MCL:
% low-order: LLF, high-order: LW; limited antidiffusion as in user's scalar code.
%
% For PN: pass fluxfun_cell handle (returns f(u) per cell), set f_cell=[]
% For MN: pass precomputed f_cell = f(u) per cell, fluxfun_cell=[]

N   = opts.N;
nM  = N+1;
Nx  = opts.Nx;
dt  = opts.dt;
dx  = opts.dx;
cfl = opts.cfl;

% Vacuum ghost state (isotropic)
ub = zeros(nM,1);
ub(1) = opts.rho_vac;

% Extend cells with ghosts
u_ext = [ub, u, ub];          % (nM x (Nx+2))

% Flux per cell
if ~isempty(f_cell)
    fb = zeros(nM,1);
    % For MN ghost: isotropic => alpha0 = log(rho/2), others 0
    if strcmpi(opts.model,'MN')
        % compute flux at ghost from isotropic alpha
        alpha_b = zeros(nM,1); alpha_b(1)=log(opts.rho_vac/2);
        % build minimal quad/basis on-the-fly is expensive; instead:
        % approximate ghost flux by zero for higher moments and use PN-like stream:
        % but for isotropic distribution, true flux is zero in all moments except moment 1?
        % Actually for isotropic psi, <mu P0 psi>=0 and <mu Pk psi> nonzero? For odd? It's 0 by symmetry.
        fb(:) = 0;
    else
        fb = flux_PN(ub, N);
    end
    f_ext = [fb, f_cell, fb]; % (nM x (Nx+2))
else
    f_ext = fluxfun_cell(u_ext); % compute in one go: accepts (nM x (Nx+2))
end

% Interface left/right states (between extended cells)
ul = u_ext(:,1:end-1);            % (nM x (Nx+1))
ur = u_ext(:,2:end);              % (nM x (Nx+1))
fl = f_ext(:,1:end-1);
fr = f_ext(:,2:end);

% Rusanov max speed bound (use 1 for slab streaming)
lambda = 1.0;                     % scalar; can be made local if you have better bounds
lam = lambda;                     % scalar

% Central flux
fCD = 0.5*(fr + fl);

% LLF (Rusanov) flux
fLF = fCD - 0.5*lam*(ur - ul);

% LW flux (user's form)
fLW = fCD - 0.5*cfl*(lam^2)*(ur - ul);

if opts.ilim == 1
    % Antidiffusive flux
    fA = fLF - fLW;

    % Local bounds per cell (componentwise) on extended cell indices
    umin_ext = u_ext;
    umax_ext = u_ext;
    if Nx >= 2
        uL = u_ext(:,1:end-2);
        uC = u_ext(:,2:end-1);
        uR = u_ext(:,3:end);
        umin_mid = min(uL, min(uC,uR));
        umax_mid = max(uL, max(uC,uR));
        umin_ext(:,2:end-1) = umin_mid;
        umax_ext(:,2:end-1) = umax_mid;
    end
    % boundary ghosts stay fixed: umin=umax=u itself

    % left/right cell bounds for each interface
    umin_L = umin_ext(:,1:end-1);
    umax_L = umax_ext(:,1:end-1);
    umin_R = umin_ext(:,2:end);
    umax_R = umax_ext(:,2:end);

    % Scaled bar states (componentwise version of your scalar formula)
    wbar = 0.5*(ur + ul)*lam - 0.5*(fr - fl);

    % Limiter (exactly your structure, componentwise)
    pos = max(0, fA);
    neg = min(0, fA);

    cap_pos = min(lam*umax_L - wbar, wbar - lam*umin_R);
    cap_neg = max(lam*umin_L - wbar, wbar - lam*umax_R);

    fA_lim = min(pos, cap_pos) + max(neg, cap_neg);

    % Flux correction
    fMCL = fLF - fA_lim;

elseif opts.ilim == -1
    fMCL = fLF;
else
    fMCL = fLW;
end

% Conservative FV update on true cells (Nx cells):
% u_i^{n+1} = u_i^n - (dt/dx) (F_{i+1/2} - F_{i-1/2})
u_new = u - (dt/dx)*( fMCL(:,2:end) - fMCL(:,1:end-1) );
end

function u = source_halfstep(u, tau, sigma_s, sigma_a, Q0)
% Exact update for isotropic scattering+absorption+isotropic source in moments.
sigma_t = sigma_s + sigma_a;

u_iso = u;
u_iso(2:end,:) = 0;

ea = exp(-sigma_a*tau);
es = exp(-sigma_s*tau);

u = ea*( es*u + (1-es)*u_iso );

if abs(Q0) > 0
    if sigma_a > 0
        u(1,:) = u(1,:) + (1-exp(-sigma_a*tau))/sigma_a * Q0;
    else
        u(1,:) = u(1,:) + tau * Q0;
    end
end
end

%% ---------------- PN flux (exact linear moment flux) ----------------
function f = flux_PN(u, N)
% u: (N+1 x Nx) moments u_l = ∫_{-1}^1 P_l(μ) ψ(μ) dμ
% flux components: f_l = ∫ μ P_l ψ dμ = (l+1)/(2l+1) u_{l+1} + l/(2l+1) u_{l-1}
if isvector(u), u = u(:); end
nM = N+1;
Nx = size(u,2);
f = zeros(nM, Nx);

% l=0
if N >= 1
    f(1,:) = u(2,:);
else
    f(1,:) = 0;
end

for l=1:N-1
    m = l+1;
    f(m,:) = ((l+1)/(2*l+1))*u(m+1,:) + (l/(2*l+1))*u(m-1,:);
end

if N >= 1
    f(N+1,:) = (N/(2*N+1))*u(N,:); % u_N = moment l=N-1
end
end

%% ---------------- MN closure (Newton for multipliers) ----------------
function [alpha_cache, ok] = multipliers_MN(u, alpha_cache, opts, Pfull, quad)
% Solve u = < b exp(alpha·b) > cellwise via damped Newton (density-scaled).
nM = size(u,1);
Nx = size(u,2);
w = quad.w_full(:);
ok = true;

for i=1:Nx
    ui = u(:,i);
    rho = ui(1);

    if rho < opts.rho_vac
        % guard -> isotropic vacuum
        ui = zeros(nM,1); ui(1)=opts.rho_vac;
        u(:,i) = ui;
        alpha_cache(:,i) = [log(opts.rho_vac/2); zeros(nM-1,1)];
        continue;
    end

    phi = ui / rho;

    beta = alpha_cache(:,i);
    beta(1) = beta(1) - log(max(rho, 1e-300));

    [beta, success] = newton_beta(phi, beta, Pfull, w, opts.newton_tol, opts.newton_maxit, opts.newton_damp);
    if ~success
        ok = false;
        % fallback isotropic (scaled)
        beta = zeros(nM,1);
        beta(1) = log(1/2);
    end

    ubeta = moment_map(beta, Pfull, w);
    rho_beta = max(ubeta(1), 1e-300);

    alpha = beta;
    alpha(1) = alpha(1) + log(rho/rho_beta);

    alpha_cache(:,i) = alpha;
end
end

function [beta, success] = newton_beta(phi, beta, Pfull, w, tol, maxit, damp)
success = true;
for it=1:maxit
    ub = moment_map(beta, Pfull, w);
    F  = ub - phi;

    if norm(F,2) < tol
        return;
    end

    H = hessian_map(beta, Pfull, w);
    d = H \ F;

    t = 1.0;
    nF = norm(F,2);
    for ls=1:20
        beta_try = beta - t*d;
        ub_try = moment_map(beta_try, Pfull, w);
        F_try = ub_try - phi;
        if isfinite(norm(F_try,2)) && norm(F_try,2) < (1 - 1e-4*t)*nF
            beta = beta_try;
            break;
        end
        t = t*damp;
    end
    if t < 1e-12
        success = false;
        return;
    end
end
success = false;
end

function u = moment_map(alpha, Pfull, w)
z = exp( (alpha.' * Pfull).');     % (nq x 1)
u = Pfull * (w .* z);              % (nM x 1)
end

function H = hessian_map(alpha, Pfull, w)
z  = exp( (alpha.' * Pfull).');    % (nq x 1)
ws = (w .* z).';                   % (1 x nq)
Bws = Pfull .* ws;                 % (nM x nq)
H = Bws * Pfull.';                 % (nM x nM)
end

function f = flux_from_alpha(alpha, Pfull, wmu_full)
% f = < μ b exp(alpha·b) >
if isvector(alpha), alpha = alpha(:); end
Nx = size(alpha,2);
psi = exp( (alpha.' * Pfull).' );     % (nq x Nx)
f = Pfull * (wmu_full .* psi);        % (nM x Nx)
end

%% ---------------- SN reference using same MCL idea (per-angle) ----------------
function rho = solve_SN_reference_MCL(ref)
Nx = ref.Nx; dx = ref.dx; dt = ref.dt; tf = ref.tf;
Nt = ceil(tf/dt); dt = tf/Nt; cfl = dt/dx;

[mu, w] = gauss_legendre(ref.Nmu, -1, 1);
mu = mu(:); w = w(:);

psi = ref.psi_vac * ones(ref.Nmu, Nx);
midL = Nx/2; midR = Nx/2 + 1;
psi(:,[midL midR]) = ref.psi_vac + 1/(2*dx);

for n=1:Nt
    psi = scatter_halfstep_SN(psi, w, dt/2, ref.sigma_s, ref.sigma_a, ref.Q0);

    % transport each angle with scalar MCL on flux f(u)=a*u
    for j=1:ref.Nmu
        a = mu(j);
        psi(j,:) = transport_scalar_MCL(psi(j,:), a, cfl, ref.ilim, ref.psi_vac);
    end

    psi = scatter_halfstep_SN(psi, w, dt/2, ref.sigma_s, ref.sigma_a, ref.Q0);
end

rho = (w.' * psi); % rho(x)=<psi>
end

function u = transport_scalar_MCL(u, a, cfl, ilim, u_vac)
% scalar conservative update on cell averages with vacuum ghosts
Nx = numel(u);
u_ext = [u_vac, u, u_vac];
ul = u_ext(1:end-1); ur = u_ext(2:end);
fl = a*ul; fr = a*ur;

lambda = abs(a);
fCD = 0.5*(fr+fl);
fLF = fCD - 0.5*lambda*(ur-ul);
fLW = fCD - 0.5*cfl*(lambda^2)*(ur-ul);

if ilim == 1
    fA = fLF - fLW;

    umin_ext = u_ext;
    umax_ext = u_ext;
    if Nx >= 2
        uL = u_ext(1:end-2); uC=u_ext(2:end-1); uR=u_ext(3:end);
        umin_ext(2:end-1) = min(uL, min(uC,uR));
        umax_ext(2:end-1) = max(uL, max(uC,uR));
    end

    umin_L = umin_ext(1:end-1); umax_L = umax_ext(1:end-1);
    umin_R = umin_ext(2:end);   umax_R = umax_ext(2:end);

    wbar = 0.5*(ur+ul)*lambda - 0.5*(fr-fl);

    pos = max(0,fA);
    neg = min(0,fA);

    cap_pos = min(lambda.*umax_L - wbar, wbar - lambda.*umin_R);
    cap_neg = max(lambda.*umin_L - wbar, wbar - lambda.*umax_R);

    fA_lim = min(pos, cap_pos) + max(neg, cap_neg);

    fMCL = fLF - fA_lim;
elseif ilim == -1
    fMCL = fLF;
else
    fMCL = fLW;
end

u = u - cfl*( fMCL(2:end) - fMCL(1:end-1) );
end

function psi = scatter_halfstep_SN(psi, w, tau, sigma_s, sigma_a, Q0)
sigma_t = sigma_s + sigma_a;

phi = (w.' * psi);
psi_iso = (phi/2);
psi_iso_mat = ones(size(psi,1),1) * psi_iso;

ea = exp(-sigma_a*tau);
es = exp(-sigma_s*tau);
psi = ea*( es*psi + (1-es)*psi_iso_mat );

if abs(Q0) > 0
    if sigma_a > 0
        add = (1-exp(-sigma_a*tau))/sigma_a * (Q0/2);
    else
        add = tau * (Q0/2);
    end
    psi = psi + add;
end
end

%% ---------------- Quadrature / Legendre utilities ----------------
function quad = make_quadrature(nq)
[mu_full, w_full] = gauss_legendre(nq, -1, 1);
quad.mu_full = mu_full(:);
quad.w_full  = w_full(:);
end

function P = legendreP_matrix(N, mu)
mu = mu(:).';
P = zeros(N+1, numel(mu));
P(1,:) = 1;
if N==0, return; end
P(2,:) = mu;
for l=1:N-1
    P(l+2,:) = ((2*l+1).*mu.*P(l+1,:) - l*P(l,:)) / (l+1);
end
end

function [x, w] = gauss_legendre(n, a, b)
i = (1:n-1)';
beta = i ./ sqrt(4*i.^2 - 1);
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x0 = diag(D);
[x0, idx] = sort(x0);
V = V(:,idx);
w0 = 2*(V(1,:).^2)';

x = (b-a)/2 * x0 + (a+b)/2;
w = (b-a)/2 * w0;
end
