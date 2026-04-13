function [u_next, step_state] = mm_step_flux_rk2(u, model, phys, grid, dt, step_cfg, cache_state)
%MM_STEP_FLUX_RK2 Second-order SSP RK2 step for the flux subsystem.
% References (Seminarquelle 2):
% - Semi-discrete flux system update: Eq. (4.13)
% - Kinetic numerical flux: Eq. (4.15)
% - Reconstruction and limiter per RK stage: Eq. (4.16)-(4.18), Lemma 4.5

if nargin < 7
    cache_state = struct();
end
tFluxTotal = tic;

if isfield(step_cfg, 'quad_flux')
    quad_flux = step_cfg.quad_flux;
elseif isfield(step_cfg, 'quad_cfg')
    quad_flux = mm_build_quadrature(model, step_cfg.quad_cfg, 'flux');
else
    quad_flux = mm_build_quadrature(model, struct(), 'flux');
end
opt_cfg = get_field_or(step_cfg, 'optimizer', struct());
rec_cfg = get_field_or(step_cfg, 'reconstruction', struct('use_characteristic', true));
lim_cfg = get_field_or(step_cfg, 'limiter', struct());
par_cfg = get_field_or(step_cfg, 'parallel', struct());
lim_type = normalize_limiter_type(get_field_or(lim_cfg, 'type', 'paper'));
if model.is_partial && model.needs_entropy && get_field_or(par_cfg, 'force_serial_partial_entropy', true)
    par_cfg.enabled = false;
end
use_char_recon = logical(get_field_or(rec_cfg, 'use_characteristic', true));
if model.is_partial && model.needs_entropy
    use_char_recon = logical(get_field_or(rec_cfg, 'use_characteristic_partial_entropy', false));
end

n_cells_top = size(u, 2);
par_cells = mm_parallel_context(par_cfg, n_cells_top, 'cells');
par_if = mm_parallel_context(par_cfg, n_cells_top + 1, 'interfaces');
use_par_cells = par_cells.use_parallel;
use_par_interfaces = par_if.use_parallel;

[tau1, s1] = flux_rhs(u, 0.0);
u1 = u + dt * tau1;

[tau2, s2] = flux_rhs(u1, dt);
u_next = u + 0.5 * dt * (tau1 + tau2);

step_state = struct();
step_state.stage1 = s1;
step_state.stage2 = s2;
step_state.cache_state = cache_state;
step_state.timing = combine_stage_timing(s1, s2, toc(tFluxTotal));

    function [rhs, rhs_state] = flux_rhs(u_stage, t_stage)
        tStage = tic;
        [nMom, nCells] = size(u_stage);

        V_cells = cell(1, nCells);
        Vinv_cells = cell(1, nCells);
        tJac = tic;
        if use_par_cells
            parfor ic = 1:nCells
                if use_char_recon
                    J = flux_jacobian_fd(u_stage(:, ic), model, quad_flux, opt_cfg);
                    [V, Vinv] = stable_characteristic_basis(J);
                    V_cells{ic} = V;
                    Vinv_cells{ic} = Vinv;
                else
                    V_cells{ic} = eye(nMom);
                    Vinv_cells{ic} = eye(nMom);
                end
            end
        else
            for ic = 1:nCells
                if use_char_recon
                    J = flux_jacobian_fd(u_stage(:, ic), model, quad_flux, opt_cfg);
                    [V, Vinv] = stable_characteristic_basis(J);
                    V_cells{ic} = V;
                    Vinv_cells{ic} = Vinv;
                else
                    V_cells{ic} = eye(nMom);
                    Vinv_cells{ic} = eye(nMom);
                end
            end
        end
        jac_state = struct();
        jac_state.V = V_cells;
        jac_state.Vinv = Vinv_cells;
        tJacobian = toc(tJac);

        tRec = tic;
        rec_cfg_apply = rec_cfg;
        rec_cfg_apply.use_characteristic = use_char_recon;
        rec_cfg_apply.parallel = par_cfg;
        [uL_rec, uR_rec, rec_state] = mm_reconstruct_characteristic(u_stage, jac_state, grid, rec_cfg_apply);
        tReconstruct = toc(tRec);

        tLimFlux = tic;
        if strcmp(lim_type, 'mcl')
            [g, lim_state] = mcl_limited_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg, use_par_interfaces);
        else
            lim_cfg_apply = lim_cfg;
            lim_cfg_apply.parallel = par_cfg;
            if get_field_or(lim_cfg, 'paper_lp_characteristic', false)
                lim_cfg_apply.characteristic_basis = jac_state;
            end
            [uL_lim, uR_lim, lim_state] = mm_apply_realizability_limiter(u_stage, uL_rec, uR_rec, model, lim_cfg_apply, quad_flux);
            g = zeros(nMom, nCells + 1);

            % Left boundary interface.
            [g(:, 1), ~] = boundary_interface_flux('left', uL_lim(:, 1), model, quad_flux, opt_cfg, phys, t_stage);

            % Interior interfaces.
            if use_par_interfaces && nCells > 1
                g_int = zeros(nMom, nCells - 1);
                parfor iif = 1:(nCells - 1)
                    [g_int(:, iif), ~] = mm_kinetic_flux(uR_lim(:, iif), uL_lim(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
                end
                g(:, 2:nCells) = g_int;
            else
                for iif = 1:(nCells - 1)
                    [g(:, iif + 1), ~] = mm_kinetic_flux(uR_lim(:, iif), uL_lim(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
                end
            end

            % Right boundary interface.
            [g(:, nCells + 1), ~] = boundary_interface_flux('right', uR_lim(:, nCells), model, quad_flux, opt_cfg, phys, t_stage);
        end
        tLimiterFlux = toc(tLimFlux);

        % Finite-volume flux difference form, Eq. (4.13).
        tRhs = tic;
        rhs = -(g(:, 2:end) - g(:, 1:end-1)) / grid.dz;
        tRhsBuild = toc(tRhs);

        rhs_state = struct();
        rhs_state.g = g;
        rhs_state.reconstruction = rec_state;
        rhs_state.limiter = lim_state;
        rhs_state.limiter_type = lim_type;
        rhs_state.timing = struct();
        rhs_state.timing.stage_total_s = toc(tStage);
        rhs_state.timing.jacobian_s = tJacobian;
        rhs_state.timing.reconstruct_s = tReconstruct;
        rhs_state.timing.limiter_and_flux_s = tLimiterFlux;
        rhs_state.timing.rhs_build_s = tRhsBuild;
    end

end

function [g, lim_state] = mcl_limited_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg, use_par_interfaces)
[nMom, nCells] = size(u_stage);
g_ho = zeros(nMom, nCells + 1);

uBLeft = boundary_moment_state('left', t_stage, model, quad_flux, phys);
uBRight = boundary_moment_state('right', t_stage, model, quad_flux, phys);

% High-order fluxes from reconstructed states.
[g_ho(:, 1), ~] = boundary_interface_flux('left', uL_rec(:, 1), model, quad_flux, opt_cfg, phys, t_stage);
if use_par_interfaces && nCells > 1
    g_ho_int = zeros(nMom, nCells - 1);
    parfor iif = 1:(nCells - 1)
        [g_ho_int(:, iif), ~] = mm_kinetic_flux(uR_rec(:, iif), uL_rec(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
    end
    g_ho(:, 2:nCells) = g_ho_int;
else
    for iif = 1:(nCells - 1)
        [g_ho(:, iif + 1), ~] = mm_kinetic_flux(uR_rec(:, iif), uL_rec(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
    end
end
[g_ho(:, nCells + 1), ~] = boundary_interface_flux('right', uR_rec(:, nCells), model, quad_flux, opt_cfg, phys, t_stage);

g_lo = zeros(nMom, nCells + 1);
g = zeros(nMom, nCells + 1);
theta = ones(1, nCells + 1);
bar_ok = true(1, nCells + 1);

lambdaMCL = get_field_or(lim_cfg, 'mcl_lambda', 1.0);
lambdaMCL = max(lambdaMCL, 1.0e-12);
maxBisect = get_field_or(lim_cfg, 'mcl_bisect_iter', 40);
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

if use_par_interfaces
    parfor iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, uBLeft, uBRight, g_ho, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect);
        g_lo(:, iif) = g_lo_i;
        g(:, iif) = g_i;
        theta(iif) = theta_i;
        bar_ok(iif) = bar_ok_i;
    end
else
    for iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, uBLeft, uBRight, g_ho, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect);
        g_lo(:, iif) = g_lo_i;
        g(:, iif) = g_i;
        theta(iif) = theta_i;
        bar_ok(iif) = bar_ok_i;
    end
end

lim_state = struct();
lim_state.theta = theta;
lim_state.mode = 'mcl';
lim_state.bar_state_feasible = bar_ok;
lim_state.g_low = g_lo;
lim_state.g_high = g_ho;

end

function [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, uBLeft, uBRight, g_ho, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect)
[~, nCells] = size(u_stage);

if iif == 1
    UL = uBLeft;
    UR = u_stage(:, 1);
elseif iif == nCells + 1
    UL = u_stage(:, nCells);
    UR = uBRight;
else
    UL = u_stage(:, iif - 1);
    UR = u_stage(:, iif);
end

[fL, ~] = mm_eval_flux_function(UL, model, quad_flux, opt_cfg, struct());
[fR, ~] = mm_eval_flux_function(UR, model, quad_flux, opt_cfg, struct());

g_lo_i = 0.5 * (fL + fR) - 0.5 * lambdaMCL * (UR - UL);
Ubar = 0.5 * (UL + UR) - (fR - fL) / (2 * lambdaMCL);
A = g_ho(:, iif) - g_lo_i;

[okBar, ~] = mm_is_realizable(Ubar, model, quad_flux, struct('epsR', epsR));
if ~okBar
    theta_i = 0.0;
    bar_ok_i = false;
    g_i = g_lo_i;
    return;
end

lo = 0.0;
hi = 1.0;
for it = 1:maxBisect
    mid = 0.5 * (lo + hi);
    Up = Ubar + (mid / lambdaMCL) * A;
    Um = Ubar - (mid / lambdaMCL) * A;

    [okP, ~] = mm_is_realizable(Up, model, quad_flux, struct('epsR', epsR));
    [okM, ~] = mm_is_realizable(Um, model, quad_flux, struct('epsR', epsR));

    if okP && okM
        lo = mid;
    else
        hi = mid;
    end
end

theta_i = lo;
bar_ok_i = true;
g_i = g_lo_i + theta_i * A;
end

function J = flux_jacobian_fd(u, model, quad_flux, opt_cfg)
n = numel(u);
J = zeros(n, n);
[f0, ~] = mm_eval_flux_function(u, model, quad_flux, opt_cfg, struct());

for j = 1:n
    e = zeros(n, 1);
    scale = max(1.0, abs(u(j)));
    epsFD = 1.0e-7 * scale;
    e(j) = epsFD;

    [fp, ~] = mm_eval_flux_function(u + e, model, quad_flux, opt_cfg, struct());
    [fm, ~] = mm_eval_flux_function(u - e, model, quad_flux, opt_cfg, struct());

    J(:, j) = (fp - fm) / (2 * epsFD);
end

if any(~isfinite(J(:)))
    J = eye(n) * max(1, norm(f0));
end
end

function [V, Vinv] = stable_characteristic_basis(J)
nLoc = size(J, 1);
V = eye(nLoc);
Vinv = eye(nLoc);

if any(~isfinite(J(:)))
    return;
end

[Vraw, D] = eig(J);
lam = diag(D);

if any(~isfinite(Vraw(:))) || any(~isfinite(lam))
    return;
end

lamScale = max(1.0, max(abs(real(lam))));
if any(abs(imag(lam)) > 1.0e-10 * lamScale)
    return;
end

[~, perm] = sort(real(lam), 'ascend');
Vtmp = real(Vraw(:, perm));

for j = 1:nLoc
    [~, imax] = max(abs(Vtmp(:, j)));
    if ~isempty(imax) && Vtmp(imax, j) < 0
        Vtmp(:, j) = -Vtmp(:, j);
    end
end

if rcond(Vtmp) < 1.0e-12 || any(~isfinite(Vtmp(:)))
    return;
end

VinvTmp = Vtmp \ eye(nLoc);
if any(~isfinite(VinvTmp(:)))
    return;
end

V = Vtmp;
Vinv = VinvTmp;
end

function [flux, uB] = boundary_interface_flux(side, uCell, model, quad_flux, opt_cfg, phys, t)
rhoB = boundary_density(side, t, phys);

if ~has_boundary_psi(side, phys)
    uB = model.b_iso * rhoB;
    if strcmp(side, 'left')
        [flux, ~] = mm_kinetic_flux(uB, uCell, model, quad_flux, struct('opt_cfg', opt_cfg));
    else
        [flux, ~] = mm_kinetic_flux(uCell, uB, model, quad_flux, struct('opt_cfg', opt_cfg));
    end
    return;
end

psiB = eval_boundary_psi(side, quad_flux.mu, t, phys, rhoB);
psiB = max(psiB(:), 0.0);
uB = quad_flux.B.' * (quad_flux.w .* psiB);

[psiCell, ~, ~] = mm_eval_ansatz(uCell, model, quad_flux, opt_cfg, struct());
psiCell = max(psiCell(:), 0.0);

if strcmp(side, 'left')
    psiInPlus = eval_boundary_psi('left', quad_flux.mu_plus, t, phys, rhoB);
    psiInPlus = max(psiInPlus(:), 0.0);
    psiOutMinus = psiCell(quad_flux.mu <= 0);
    flux = quad_flux.B_plus.' * (quad_flux.w_plus .* quad_flux.mu_plus .* psiInPlus) + ...
           quad_flux.B_minus.' * (quad_flux.w_minus .* quad_flux.mu_minus .* psiOutMinus);
else
    psiOutPlus = psiCell(quad_flux.mu >= 0);
    psiInMinus = eval_boundary_psi('right', quad_flux.mu_minus, t, phys, rhoB);
    psiInMinus = max(psiInMinus(:), 0.0);
    flux = quad_flux.B_plus.' * (quad_flux.w_plus .* quad_flux.mu_plus .* psiOutPlus) + ...
           quad_flux.B_minus.' * (quad_flux.w_minus .* quad_flux.mu_minus .* psiInMinus);
end
end

function uB = boundary_moment_state(side, t, model, quad_flux, phys)
rhoB = boundary_density(side, t, phys);
if ~has_boundary_psi(side, phys)
    uB = model.b_iso * rhoB;
    return;
end

psiB = eval_boundary_psi(side, quad_flux.mu, t, phys, rhoB);
psiB = max(psiB(:), 0.0);
uB = quad_flux.B.' * (quad_flux.w .* psiB);
end

function tf = has_boundary_psi(side, phys)
tf = false;
if ~isfield(phys, 'boundary_psi')
    return;
end

bpsi = phys.boundary_psi;
if isa(bpsi, 'function_handle')
    tf = true;
elseif isstruct(bpsi)
    if strcmp(side, 'left')
        tf = isfield(bpsi, 'left');
    else
        tf = isfield(bpsi, 'right');
    end
end
end

function psi = eval_boundary_psi(side, mu, t, phys, default_rho)
psi = default_rho * ones(size(mu));
if ~isfield(phys, 'boundary_psi')
    return;
end

bpsi = phys.boundary_psi;
if isa(bpsi, 'function_handle')
    psi = call_boundary_psi(bpsi, side, mu, t);
elseif isstruct(bpsi)
    if strcmp(side, 'left') && isfield(bpsi, 'left')
        psi = eval_boundary_psi_item(bpsi.left, side, mu, t);
    elseif strcmp(side, 'right') && isfield(bpsi, 'right')
        psi = eval_boundary_psi_item(bpsi.right, side, mu, t);
    end
end

if isscalar(psi)
    psi = psi * ones(size(mu));
end
psi = reshape(psi, size(mu));
end

function psi = eval_boundary_psi_item(item, side, mu, t)
if isa(item, 'function_handle')
    psi = call_boundary_psi(item, side, mu, t);
else
    psi = item;
end
end

function psi = call_boundary_psi(fun, side, mu, t)
try
    psi = fun(side, mu, t);
    return;
catch
end
try
    psi = fun(mu, t);
    return;
catch
end
try
    psi = fun(mu);
    return;
catch
end
try
    psi = fun(t, mu);
    return;
catch
end
try
    psi = fun(t);
catch
    error('Unsupported boundary_psi function signature.');
end
end

function rho = boundary_density(side, t, phys)
if isfield(phys, 'boundary')
    b = phys.boundary;
    if isa(b, 'function_handle')
        rho = b(side, t);
        return;
    elseif isstruct(b)
        if strcmp(side, 'left') && isfield(b, 'left')
            rho = eval_boundary_item(b.left, t);
            return;
        elseif strcmp(side, 'right') && isfield(b, 'right')
            rho = eval_boundary_item(b.right, t);
            return;
        end
    end
end

if isfield(phys, 'psi_vac_density')
    rho = phys.psi_vac_density;
else
    rho = 0.0;
end
end

function rho = eval_boundary_item(item, t)
if isa(item, 'function_handle')
    rho = item(t);
else
    rho = item;
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function out = normalize_limiter_type(in)
if isstring(in)
    out = lower(char(in));
elseif ischar(in)
    out = lower(strtrim(in));
else
    out = 'paper';
end
if isempty(out)
    out = 'paper';
end
end

function timing = combine_stage_timing(s1, s2, totalTime)
timing = struct();
timing.total_s = totalTime;
timing.stage1_total_s = get_timing_field(s1, 'stage_total_s');
timing.stage2_total_s = get_timing_field(s2, 'stage_total_s');
timing.jacobian_s = get_timing_field(s1, 'jacobian_s') + get_timing_field(s2, 'jacobian_s');
timing.reconstruct_s = get_timing_field(s1, 'reconstruct_s') + get_timing_field(s2, 'reconstruct_s');
timing.limiter_and_flux_s = get_timing_field(s1, 'limiter_and_flux_s') + get_timing_field(s2, 'limiter_and_flux_s');
timing.rhs_build_s = get_timing_field(s1, 'rhs_build_s') + get_timing_field(s2, 'rhs_build_s');
end

function v = get_timing_field(stageState, name)
v = 0.0;
if isfield(stageState, 'timing') && isfield(stageState.timing, name)
    v = stageState.timing.(name);
end
end
