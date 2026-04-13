function [state_next, step_state] = mm_step_holo_m1pn_flux_rk2(state, models, phys, grid, dt, step_cfg, cache_state)
%MM_STEP_HOLO_M1PN_FLUX_RK2 Coupled RK2 flux step for HOLO M1+PN scheme.
% Base building blocks reuse Seminarquelle 2 equations:
% - Flux difference update Eq. (4.13)
% - Kinetic interface flux Eq. (4.15)
% - Reconstruction/limiting Eq. (4.16)-(4.18)

if nargin < 7 || isempty(cache_state)
    cache_state = struct();
end
tFluxTotal = tic;

model_pn = models.pn;
model_m1 = models.m1;

u0_pn = state.u_pn;
if isfield(state, 'u_m1') && ~isempty(state.u_m1)
    u0_m1 = state.u_m1;
else
    u0_m1 = u0_pn(1:2, :);
end

if size(u0_m1, 1) ~= 2
    error('M1 state must have exactly 2 moments.');
end
if size(u0_pn, 2) ~= size(u0_m1, 2)
    error('PN and M1 states must use the same number of cells.');
end

if isfield(step_cfg, 'quad_flux_pn')
    quad_pn = step_cfg.quad_flux_pn;
elseif isfield(step_cfg, 'quad_cfg')
    quad_pn = mm_build_quadrature(model_pn, step_cfg.quad_cfg, 'flux');
else
    quad_pn = mm_build_quadrature(model_pn, struct(), 'flux');
end

if isfield(step_cfg, 'quad_flux_m1')
    quad_m1_flux = step_cfg.quad_flux_m1;
elseif isfield(step_cfg, 'quad_cfg')
    quad_m1_flux = mm_build_quadrature(model_m1, step_cfg.quad_cfg, 'flux');
else
    quad_m1_flux = mm_build_quadrature(model_m1, struct(), 'flux');
end

if isfield(step_cfg, 'quad_lp_m1')
    quad_m1_lp = step_cfg.quad_lp_m1;
elseif isfield(step_cfg, 'quad_cfg')
    quad_m1_lp = mm_build_quadrature(model_m1, step_cfg.quad_cfg, 'lp');
else
    quad_m1_lp = mm_build_quadrature(model_m1, struct(), 'lp');
end

opt_cfg = get_field_or(step_cfg, 'optimizer', struct());
rec_cfg = get_field_or(step_cfg, 'reconstruction', struct('use_characteristic', true));
lim_cfg = get_field_or(step_cfg, 'limiter', struct());
holo_cfg = get_field_or(step_cfg, 'holo', struct());
par_cfg = get_field_or(step_cfg, 'parallel', struct());

n_cells_top = size(u0_pn, 2);
par_cells = mm_parallel_context(par_cfg, n_cells_top, 'cells');
par_if = mm_parallel_context(par_cfg, n_cells_top + 1, 'interfaces');
use_par_cells = par_cells.use_parallel;
use_par_interfaces = par_if.use_parallel;

use_mcl_pn = logical(get_field_or(holo_cfg, 'use_mcl_pn', true));
use_mcl_m1 = logical(get_field_or(holo_cfg, 'use_mcl_m1', true));

sync_cfg = get_field_or(step_cfg, 'sync', struct());
sync_cfg.weight = get_field_or(sync_cfg, 'weight', get_field_or(holo_cfg, 'sync_weight', 'identity'));
sync_each_stage = logical(get_field_or(sync_cfg, 'each_stage', get_field_or(holo_cfg, 'sync_each_stage', true)));

[k1_pn, k1_m1, s1] = coupled_rhs(u0_pn, u0_m1, 0.0);
u1_pn_ho = u0_pn + dt * k1_pn;
u1_m1 = u0_m1 + dt * k1_m1;

if sync_each_stage
    [u1_pn, sync_stage1] = mm_sync_pn_to_m1(u1_pn_ho, u1_m1, sync_cfg);
    u1_m1 = u1_pn(1:2, :);
else
    u1_pn = u1_pn_ho;
    sync_stage1 = empty_sync_state(u1_pn, u1_m1);
end

[k2_pn, k2_m1, s2] = coupled_rhs(u1_pn, u1_m1, dt);
u2_pn_ho = u0_pn + 0.5 * dt * (k1_pn + k2_pn);
u2_m1 = u0_m1 + 0.5 * dt * (k1_m1 + k2_m1);

[u2_pn, sync_final] = mm_sync_pn_to_m1(u2_pn_ho, u2_m1, sync_cfg);
u2_m1 = u2_pn(1:2, :);

state_next = struct();
state_next.u_pn = u2_pn;
state_next.u_m1 = u2_m1;
state_next.cache_pn = get_field_or(cache_state, 'cache_pn', struct());
state_next.cache_m1 = get_field_or(cache_state, 'cache_m1', struct());
state_next.diag = struct('sync_stage1', sync_stage1, 'sync_final', sync_final);

step_state = struct();
step_state.stage1 = s1;
step_state.stage2 = s2;
step_state.sync_stage1 = sync_stage1;
step_state.sync_final = sync_final;
step_state.cache_state = cache_state;
step_state.diag = struct();
step_state.diag.theta_pn_active_stage1 = count_active(get_field_or(s1.pn_limiter, 'theta', []));
step_state.diag.theta_pn_active_stage2 = count_active(get_field_or(s2.pn_limiter, 'theta', []));
step_state.diag.theta_m1_active_stage1 = count_active(get_field_or(s1.m1_limiter, 'theta', []));
step_state.diag.theta_m1_active_stage2 = count_active(get_field_or(s2.m1_limiter, 'theta', []));
step_state.diag.sync_corr_stage1 = sync_stage1.correction_norm_fro;
step_state.diag.sync_corr_final = sync_final.correction_norm_fro;
step_state.timing = combine_holo_stage_timing(s1, s2, sync_stage1, sync_final, toc(tFluxTotal));

    function [rhs_pn, rhs_m1, rhs_state] = coupled_rhs(u_pn_stage, u_m1_stage, t_stage)
        tStage = tic;
        [nMomPn, nCells] = size(u_pn_stage);

        V_cells = cell(1, nCells);
        Vinv_cells = cell(1, nCells);
        tJac = tic;
        if use_par_cells
            parfor ic = 1:nCells
                if rec_cfg.use_characteristic
                    J = flux_jacobian_fd(u_pn_stage(:, ic), model_pn, quad_pn, opt_cfg);
                    [V, D] = eig(J);
                    if any(~isfinite(V(:))) || any(~isfinite(diag(D)))
                        V = eye(nMomPn);
                    end
                    V_cells{ic} = V;
                    Vinv_cells{ic} = pinv(V);
                else
                    V_cells{ic} = eye(nMomPn);
                    Vinv_cells{ic} = eye(nMomPn);
                end
            end
        else
            for ic = 1:nCells
                if rec_cfg.use_characteristic
                    J = flux_jacobian_fd(u_pn_stage(:, ic), model_pn, quad_pn, opt_cfg);
                    [V, D] = eig(J);
                    if any(~isfinite(V(:))) || any(~isfinite(diag(D)))
                        V = eye(nMomPn);
                    end
                    V_cells{ic} = V;
                    Vinv_cells{ic} = pinv(V);
                else
                    V_cells{ic} = eye(nMomPn);
                    Vinv_cells{ic} = eye(nMomPn);
                end
            end
        end
        jac_state = struct();
        jac_state.V = V_cells;
        jac_state.Vinv = Vinv_cells;
        tJacobian = toc(tJac);

        tRec = tic;
        rec_cfg_apply = rec_cfg;
        rec_cfg_apply.parallel = par_cfg;
        [uL_rec, uR_rec, rec_state] = mm_reconstruct_characteristic(u_pn_stage, jac_state, grid, rec_cfg_apply);
        tReconstruct = toc(tRec);

        tPn = tic;
        if use_mcl_pn
            [g_pn, pn_lim_state] = mcl_interface_flux(u_pn_stage, uL_rec, uR_rec, model_pn, phys, t_stage, quad_pn, opt_cfg, lim_cfg, use_par_interfaces);
        else
            lim_cfg_apply = lim_cfg;
            lim_cfg_apply.parallel = par_cfg;
            [uL_lim, uR_lim, pn_lim_state] = mm_apply_realizability_limiter(u_pn_stage, uL_rec, uR_rec, model_pn, lim_cfg_apply, quad_pn);
            g_pn = interface_flux_from_reconstructed(uL_lim, uR_lim, model_pn, phys, t_stage, quad_pn, opt_cfg, use_par_interfaces);
            if ~isfield(pn_lim_state, 'mode')
                pn_lim_state.mode = 'paper';
            end
        end
        tPnFlux = toc(tPn);

        tRhsPn = tic;
        rhs_pn = -(g_pn(:, 2:end) - g_pn(:, 1:end-1)) / grid.dz;
        tRhsPnBuild = toc(tRhsPn);

        g_ho_m1 = g_pn(1:2, :);
        tM1 = tic;
        [g_m1, m1_lim_state] = m1_holo_mcl_flux(u_m1_stage, g_ho_m1, model_m1, phys, t_stage, quad_m1_flux, quad_m1_lp, opt_cfg, lim_cfg, use_mcl_m1, use_par_interfaces);
        tM1Flux = toc(tM1);
        tRhsM1 = tic;
        rhs_m1 = -(g_m1(:, 2:end) - g_m1(:, 1:end-1)) / grid.dz;
        tRhsM1Build = toc(tRhsM1);

        rhs_state = struct();
        rhs_state.g_pn = g_pn;
        rhs_state.g_m1 = g_m1;
        rhs_state.reconstruction = rec_state;
        rhs_state.pn_limiter = pn_lim_state;
        rhs_state.m1_limiter = m1_lim_state;
        rhs_state.timing = struct();
        rhs_state.timing.stage_total_s = toc(tStage);
        rhs_state.timing.jacobian_s = tJacobian;
        rhs_state.timing.reconstruct_s = tReconstruct;
        rhs_state.timing.pn_limiter_and_flux_s = tPnFlux;
        rhs_state.timing.pn_rhs_build_s = tRhsPnBuild;
        rhs_state.timing.m1_limiter_and_flux_s = tM1Flux;
        rhs_state.timing.m1_rhs_build_s = tRhsM1Build;
    end

end

function [g, lim_state] = m1_holo_mcl_flux(u_stage, g_ho, model_m1, phys, t_stage, quad_flux, quad_lp, opt_cfg, lim_cfg, use_mcl, use_par_interfaces)
[~, nCells] = size(u_stage);

lambdaMCL = get_field_or(lim_cfg, 'mcl_lambda', 1.0);
lambdaMCL = max(lambdaMCL, 1.0e-12);
maxBisect = get_field_or(lim_cfg, 'mcl_bisect_iter', 40);
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

g_lo = zeros(2, nCells + 1);
g = zeros(2, nCells + 1);
theta = ones(1, nCells + 1);
bar_ok = true(1, nCells + 1);

uBLeft = model_m1.b_iso * boundary_density('left', t_stage, phys);
uBRight = model_m1.b_iso * boundary_density('right', t_stage, phys);

if use_par_interfaces
    parfor iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = m1_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model_m1, quad_flux, quad_lp, opt_cfg, lambdaMCL, epsR, maxBisect, use_mcl);
        g_lo(:, iif) = g_lo_i;
        g(:, iif) = g_i;
        theta(iif) = theta_i;
        bar_ok(iif) = bar_ok_i;
    end
else
    for iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = m1_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model_m1, quad_flux, quad_lp, opt_cfg, lambdaMCL, epsR, maxBisect, use_mcl);
        g_lo(:, iif) = g_lo_i;
        g(:, iif) = g_i;
        theta(iif) = theta_i;
        bar_ok(iif) = bar_ok_i;
    end
end

lim_state = struct();
lim_state.mode = 'mcl';
lim_state.theta = theta;
lim_state.bar_state_feasible = bar_ok;
lim_state.g_low = g_lo;
lim_state.g_high = g_ho;

end

function [g_lo_i, g_i, theta_i, bar_ok_i] = m1_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model_m1, quad_flux, quad_lp, opt_cfg, lambdaMCL, epsR, maxBisect, use_mcl)
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

[fL, ~] = mm_eval_flux_function(UL, model_m1, quad_flux, opt_cfg, struct());
[fR, ~] = mm_eval_flux_function(UR, model_m1, quad_flux, opt_cfg, struct());
g_lo_i = 0.5 * (fL + fR) - 0.5 * lambdaMCL * (UR - UL);

if ~use_mcl
    g_i = g_ho(:, iif);
    theta_i = 1.0;
    bar_ok_i = true;
    return;
end

Ubar = 0.5 * (UL + UR) - (fR - fL) / (2 * lambdaMCL);
A = g_ho(:, iif) - g_lo_i;

[okBar, ~] = mm_is_realizable(Ubar, model_m1, quad_lp, struct('epsR', epsR));
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

    [okP, ~] = mm_is_realizable(Up, model_m1, quad_lp, struct('epsR', epsR));
    [okM, ~] = mm_is_realizable(Um, model_m1, quad_lp, struct('epsR', epsR));

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

function [g, lim_state] = mcl_interface_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg, use_par_interfaces)
[nMom, nCells] = size(u_stage);
g_ho = interface_flux_from_reconstructed(uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, use_par_interfaces);

g_lo = zeros(nMom, nCells + 1);
g = zeros(nMom, nCells + 1);
theta = ones(1, nCells + 1);
bar_ok = true(1, nCells + 1);

lambdaMCL = get_field_or(lim_cfg, 'mcl_lambda', 1.0);
lambdaMCL = max(lambdaMCL, 1.0e-12);
maxBisect = get_field_or(lim_cfg, 'mcl_bisect_iter', 40);
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

uBLeft = model.b_iso * boundary_density('left', t_stage, phys);
uBRight = model.b_iso * boundary_density('right', t_stage, phys);

if use_par_interfaces
    parfor iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect);
        g_lo(:, iif) = g_lo_i;
        g(:, iif) = g_i;
        theta(iif) = theta_i;
        bar_ok(iif) = bar_ok_i;
    end
else
    for iif = 1:(nCells + 1)
        [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect);
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

function [g_lo_i, g_i, theta_i, bar_ok_i] = mcl_interface_state(iif, u_stage, g_ho, uBLeft, uBRight, model, quad_flux, opt_cfg, lambdaMCL, epsR, maxBisect)
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

function g = interface_flux_from_reconstructed(uL, uR, model, phys, t_stage, quad_flux, opt_cfg, use_par_interfaces)
[nMom, nCells] = size(uL);
g = zeros(nMom, nCells + 1);

uBLeft = model.b_iso * boundary_density('left', t_stage, phys);
uBRight = model.b_iso * boundary_density('right', t_stage, phys);

[g(:, 1), ~] = mm_kinetic_flux(uBLeft, uL(:, 1), model, quad_flux, struct('opt_cfg', opt_cfg));
if use_par_interfaces && nCells > 1
    g_int = zeros(nMom, nCells - 1);
    parfor iif = 1:(nCells - 1)
        [g_int(:, iif), ~] = mm_kinetic_flux(uR(:, iif), uL(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
    end
    g(:, 2:nCells) = g_int;
else
    for iif = 1:(nCells - 1)
        [g(:, iif + 1), ~] = mm_kinetic_flux(uR(:, iif), uL(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
    end
end
[g(:, nCells + 1), ~] = mm_kinetic_flux(uR(:, nCells), uBRight, model, quad_flux, struct('opt_cfg', opt_cfg));

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

function n = count_active(theta)
if isempty(theta)
    n = 0;
else
    n = sum(theta < (1.0 - 1.0e-12));
end
end

function sync_state = empty_sync_state(u_pn, u_m1)
sync_state = struct();
sync_state.mode = 'none';
sync_state.correction_norm_fro = 0.0;
sync_state.correction_norm_max = 0.0;
sync_state.max_constraint_violation = max(abs(u_pn(1:2, :) - u_m1), [], 'all');
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function timing = combine_holo_stage_timing(s1, s2, syncStage1, syncFinal, totalTime)
timing = struct();
timing.total_s = totalTime;
timing.stage1_total_s = get_timing_field(s1, 'stage_total_s');
timing.stage2_total_s = get_timing_field(s2, 'stage_total_s');
timing.jacobian_s = get_timing_field(s1, 'jacobian_s') + get_timing_field(s2, 'jacobian_s');
timing.reconstruct_s = get_timing_field(s1, 'reconstruct_s') + get_timing_field(s2, 'reconstruct_s');
timing.pn_limiter_and_flux_s = get_timing_field(s1, 'pn_limiter_and_flux_s') + get_timing_field(s2, 'pn_limiter_and_flux_s');
timing.pn_rhs_build_s = get_timing_field(s1, 'pn_rhs_build_s') + get_timing_field(s2, 'pn_rhs_build_s');
timing.m1_limiter_and_flux_s = get_timing_field(s1, 'm1_limiter_and_flux_s') + get_timing_field(s2, 'm1_limiter_and_flux_s');
timing.m1_rhs_build_s = get_timing_field(s1, 'm1_rhs_build_s') + get_timing_field(s2, 'm1_rhs_build_s');
timing.sync_stage1_s = get_field_or(syncStage1, 'elapsed_s', 0.0);
timing.sync_final_s = get_field_or(syncFinal, 'elapsed_s', 0.0);
timing.sync_total_s = timing.sync_stage1_s + timing.sync_final_s;
end

function v = get_timing_field(stageState, name)
v = 0.0;
if isfield(stageState, 'timing') && isfield(stageState.timing, name)
    v = stageState.timing.(name);
end
end
