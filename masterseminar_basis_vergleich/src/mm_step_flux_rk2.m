function [u_next, step_state] = mm_step_flux_rk2(u, model, phys, grid, dt, step_cfg, cache_state)
%MM_STEP_FLUX_RK2 Second-order SSP RK2 step for the flux subsystem.
% References (Seminarquelle 2):
% - Semi-discrete flux system update: Eq. (4.13)
% - Kinetic numerical flux: Eq. (4.15)
% - Reconstruction and limiter per RK stage: Eq. (4.16)-(4.18), Lemma 4.5

if nargin < 7 || isempty(cache_state)
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

incoming_cache = extract_flux_cache_state(cache_state);

[tau1, s1, stage1_cache] = flux_rhs(u, 0.0, get_stage_warm_cache(incoming_cache, struct()));
u1 = u + dt * tau1;

[tau2, s2, stage2_cache] = flux_rhs(u1, dt, get_stage_warm_cache(incoming_cache, stage1_cache));
u_next = u + 0.5 * dt * (tau1 + tau2);

step_state = struct();
step_state.stage1 = s1;
step_state.stage2 = s2;
step_state.cache_state = build_next_flux_cache(incoming_cache, stage1_cache, stage2_cache);
step_state.timing = combine_stage_timing(s1, s2, toc(tFluxTotal));

    function [rhs, rhs_state, stage_cache] = flux_rhs(u_stage, t_stage, warm_stage_cache)
        tStage = tic;
        [nMom, nCells] = size(u_stage);

        tCache = tic;
        [stage_cache, jac_state, cache_prep] = prepare_stage_cell_cache( ...
            u_stage, model, quad_flux, opt_cfg, warm_stage_cache, use_char_recon, use_par_cells);
        tStageCache = toc(tCache);

        tRec = tic;
        rec_cfg_apply = rec_cfg;
        rec_cfg_apply.use_characteristic = use_char_recon;
        rec_cfg_apply.parallel = par_cfg;
        [uL_rec, uR_rec, rec_state] = mm_reconstruct_characteristic(u_stage, jac_state, grid, rec_cfg_apply);
        tReconstruct = toc(tRec);

        boundary_cache = prepare_boundary_cache(phys, t_stage, model, quad_flux);

        tLimFlux = tic;
        if strcmp(lim_type, 'mcl')
            [g, lim_state, flux_diag] = mcl_limited_flux( ...
                u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, ...
                lim_cfg, use_par_interfaces, stage_cache, boundary_cache);
        else
            lim_cfg_apply = lim_cfg;
            lim_cfg_apply.parallel = par_cfg;
            if get_field_or(lim_cfg, 'paper_lp_characteristic', false)
                lim_cfg_apply.characteristic_basis = jac_state;
            end
            [uL_lim, uR_lim, lim_state] = mm_apply_realizability_limiter(u_stage, uL_rec, uR_rec, model, lim_cfg_apply, quad_flux);

            flux_diag = struct('entropy_reconstructed_solves_s', 0.0, ...
                'high_order_flux_s', 0.0, 'low_order_flux_s', 0.0, ...
                'mcl_bisection_s', 0.0);
            entropy_rec_state = struct('available', false, 'failed_cells', false(1, nCells), ...
                'n_failed_cells', 0, 'left_alpha', [], 'right_alpha', []);

            if model.needs_entropy
                tRecEntropy = tic;
                [uL_lim, uR_lim, entropy_rec_state] = ...
                    prepare_entropy_reconstructed_states(u_stage, uL_lim, uR_lim, model, quad_flux, opt_cfg, stage_cache);
                flux_diag.entropy_reconstructed_solves_s = toc(tRecEntropy);
            end

            tHigh = tic;
            g = compute_interface_fluxes(uL_lim, uR_lim, model, phys, t_stage, quad_flux, ...
                opt_cfg, use_par_interfaces, entropy_rec_state, boundary_cache);
            flux_diag.high_order_flux_s = toc(tHigh);

            g = enforce_reflection_flux_symmetry(g, u_stage, uL_lim, uR_lim, model);
            lim_state.reconstructed_entropy_fallback_cells = entropy_rec_state.failed_cells;
            lim_state.reconstructed_entropy_fallback_count = entropy_rec_state.n_failed_cells;
        end
        tLimiterFlux = toc(tLimFlux);

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
        rhs_state.timing.stage_cache_preparation_s = tStageCache;
        rhs_state.timing.entropy_cell_solves_s = cache_prep.entropy_cell_solves_s;
        rhs_state.timing.jacobian_s = cache_prep.characteristic_basis_s;
        rhs_state.timing.characteristic_basis_s = cache_prep.characteristic_basis_s;
        rhs_state.timing.reconstruct_s = tReconstruct;
        rhs_state.timing.entropy_reconstructed_solves_s = flux_diag.entropy_reconstructed_solves_s;
        rhs_state.timing.high_order_flux_s = flux_diag.high_order_flux_s;
        rhs_state.timing.low_order_flux_s = flux_diag.low_order_flux_s;
        rhs_state.timing.mcl_bisection_s = flux_diag.mcl_bisection_s;
        rhs_state.timing.limiter_and_flux_s = tLimiterFlux;
        rhs_state.timing.rhs_build_s = tRhsBuild;
    end

end

function g = enforce_reflection_flux_symmetry(g, u_stage, uL_lim, uR_lim, model)
if ~isfield(model, 'reflection_matrix') || isempty(model.reflection_matrix)
    return;
end

R = model.reflection_matrix;
nCells = size(uL_lim, 2);

% If the whole stage state is mirrored, project the interface fluxes back
% onto the exact reflection-symmetric subspace to remove solver drift.
if stage_is_mirrored(u_stage, R)
    for ig = 2:nCells
        jg = nCells + 2 - ig;
        if ig > jg
            continue;
        end

        gi = g(:, ig);
        if ig == jg
            g(:, ig) = 0.5 * (gi - R * gi);
        else
            gj = g(:, jg);
            gi_sym = 0.5 * (gi - R * gj);
            g(:, ig) = gi_sym;
            g(:, jg) = -R * gi_sym;
        end
    end
    return;
end

for iif = 1:(nCells - 1)
    jif = nCells - iif;
    if iif > jif
        continue;
    end

    if ~interfaces_are_mirrored(uL_lim, uR_lim, R, iif, jif)
        continue;
    end

    gi = g(:, iif + 1);
    if iif == jif
        g(:, iif + 1) = 0.5 * (gi - R * gi);
    else
        gj = g(:, jif + 1);
        gi_sym = 0.5 * (gi - R * gj);
        g(:, iif + 1) = gi_sym;
        g(:, jif + 1) = -R * gi_sym;
    end
end
end

function tf = stage_is_mirrored(u_stage, R)
nCells = size(u_stage, 2);
maxErr = 0.0;
scale = max(1.0, norm(u_stage, Inf));

for ic = 1:nCells
    jc = nCells + 1 - ic;
    maxErr = max(maxErr, norm(u_stage(:, ic) - R * u_stage(:, jc), Inf));
    if maxErr > 1.0e-8 * scale
        tf = false;
        return;
    end
end

tf = true;
end

function tf = interfaces_are_mirrored(uL_lim, uR_lim, R, iif, jif)
leftErr = norm(uR_lim(:, iif) - R * uL_lim(:, jif + 1), Inf);
rightErr = norm(uL_lim(:, iif + 1) - R * uR_lim(:, jif), Inf);

scale = max(1.0, max([ ...
    norm(uR_lim(:, iif), Inf), ...
    norm(uL_lim(:, iif + 1), Inf), ...
    norm(uR_lim(:, jif), Inf), ...
    norm(uL_lim(:, jif + 1), Inf)]));

tf = leftErr <= 1.0e-8 * scale && rightErr <= 1.0e-8 * scale;
end

function [g, lim_state, flux_timing] = mcl_limited_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg, use_par_interfaces, stage_cache, boundary_cache)
[nMom, nCells] = size(u_stage);
g_ho = zeros(nMom, nCells + 1);

flux_timing = struct('entropy_reconstructed_solves_s', 0.0, ...
    'high_order_flux_s', 0.0, 'low_order_flux_s', 0.0, 'mcl_bisection_s', 0.0);

uL_ho = uL_rec;
uR_ho = uR_rec;
entropy_rec_state = struct('available', false, 'failed_cells', false(1, nCells), ...
    'n_failed_cells', 0, 'left_alpha', [], 'right_alpha', []);
if model.needs_entropy
    tRecEntropy = tic;
    [uL_ho, uR_ho, entropy_rec_state] = ...
        prepare_entropy_reconstructed_states(u_stage, uL_ho, uR_ho, model, quad_flux, opt_cfg, stage_cache);
    flux_timing.entropy_reconstructed_solves_s = toc(tRecEntropy);
end

tHigh = tic;
g_ho = compute_interface_fluxes(uL_ho, uR_ho, model, phys, t_stage, quad_flux, ...
    opt_cfg, use_par_interfaces, entropy_rec_state, boundary_cache);
flux_timing.high_order_flux_s = toc(tHigh);

uBLeft = boundary_cache.left.uB;
uBRight = boundary_cache.right.uB;
fBLeft = boundary_cache.left.fB;
fBRight = boundary_cache.right.fB;
f_cell = stage_cache.flux_cell;

g_lo = zeros(nMom, nCells + 1);
g = zeros(nMom, nCells + 1);
theta = ones(1, nCells + 1);
bar_ok = true(1, nCells + 1);
Ubar = zeros(nMom, nCells + 1);
A = zeros(nMom, nCells + 1);

lambdaMCL = get_field_or(lim_cfg, 'mcl_lambda', 1.0);
lambdaMCL = max(lambdaMCL, 1.0e-12);
maxBisect = get_field_or(lim_cfg, 'mcl_bisect_iter', 20);
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

tLow = tic;
if use_par_interfaces
    parfor iif = 1:(nCells + 1)
        [g_lo_i, Ubar_i, A_i, bar_ok_i] = build_mcl_interface_data( ...
            iif, u_stage, uBLeft, uBRight, f_cell, fBLeft, fBRight, g_ho(:, iif), ...
            model, quad_flux, lambdaMCL, epsR);
        g_lo(:, iif) = g_lo_i;
        Ubar(:, iif) = Ubar_i;
        A(:, iif) = A_i;
        bar_ok(iif) = bar_ok_i;
    end
else
    for iif = 1:(nCells + 1)
        [g_lo_i, Ubar_i, A_i, bar_ok_i] = build_mcl_interface_data( ...
            iif, u_stage, uBLeft, uBRight, f_cell, fBLeft, fBRight, g_ho(:, iif), ...
            model, quad_flux, lambdaMCL, epsR);
        g_lo(:, iif) = g_lo_i;
        Ubar(:, iif) = Ubar_i;
        A(:, iif) = A_i;
        bar_ok(iif) = bar_ok_i;
    end
end
flux_timing.low_order_flux_s = toc(tLow);

tBisect = tic;
if use_par_interfaces
    parfor iif = 1:(nCells + 1)
        [g_i, theta_i] = finalize_mcl_interface(Ubar(:, iif), A(:, iif), g_lo(:, iif), ...
            bar_ok(iif), model, quad_flux, lambdaMCL, epsR, maxBisect);
        g(:, iif) = g_i;
        theta(iif) = theta_i;
    end
else
    for iif = 1:(nCells + 1)
        [g_i, theta_i] = finalize_mcl_interface(Ubar(:, iif), A(:, iif), g_lo(:, iif), ...
            bar_ok(iif), model, quad_flux, lambdaMCL, epsR, maxBisect);
        g(:, iif) = g_i;
        theta(iif) = theta_i;
    end
end
flux_timing.mcl_bisection_s = toc(tBisect);

lim_state = struct();
lim_state.theta = theta;
lim_state.mode = 'mcl';
lim_state.bar_state_feasible = bar_ok;
lim_state.g_low = g_lo;
lim_state.g_high = g_ho;
lim_state.reconstructed_entropy_fallback_cells = entropy_rec_state.failed_cells;
lim_state.reconstructed_entropy_fallback_count = entropy_rec_state.n_failed_cells;
end

function [g_lo_i, Ubar_i, A_i, bar_ok_i] = build_mcl_interface_data(iif, u_stage, uBLeft, uBRight, f_cell, fBLeft, fBRight, g_ho_i, model, quad_flux, lambdaMCL, epsR)
[~, nCells] = size(u_stage);

if iif == 1
    UL = uBLeft;
    UR = u_stage(:, 1);
    fL = fBLeft;
    fR = f_cell(:, 1);
elseif iif == nCells + 1
    UL = u_stage(:, nCells);
    UR = uBRight;
    fL = f_cell(:, nCells);
    fR = fBRight;
else
    UL = u_stage(:, iif - 1);
    UR = u_stage(:, iif);
    fL = f_cell(:, iif - 1);
    fR = f_cell(:, iif);
end

g_lo_i = 0.5 * (fL + fR) - 0.5 * lambdaMCL * (UR - UL);
Ubar_i = 0.5 * (UL + UR) - (fR - fL) / (2 * lambdaMCL);
A_i = g_ho_i - g_lo_i;

[bar_ok_i, ~] = mm_is_realizable(Ubar_i, model, quad_flux, struct('epsR', epsR));
end

function [g_i, theta_i] = finalize_mcl_interface(Ubar_i, A_i, g_lo_i, bar_ok_i, model, quad_flux, lambdaMCL, epsR, maxBisect)
if ~bar_ok_i
    theta_i = 0.0;
    g_i = g_lo_i;
    return;
end

lo = 0.0;
hi = 1.0;
for it = 1:maxBisect
    mid = 0.5 * (lo + hi);
    Up = Ubar_i + (mid / lambdaMCL) * A_i;
    Um = Ubar_i - (mid / lambdaMCL) * A_i;

    [okP, ~] = mm_is_realizable(Up, model, quad_flux, struct('epsR', epsR));
    [okM, ~] = mm_is_realizable(Um, model, quad_flux, struct('epsR', epsR));

    if okP && okM
        lo = mid;
    else
        hi = mid;
    end
end

theta_i = lo;
g_i = g_lo_i + theta_i * A_i;
end

function g = compute_interface_fluxes(uL, uR, model, phys, t_stage, quad_flux, opt_cfg, use_par_interfaces, entropy_state, boundary_cache)
[nMom, nCells] = size(uL);
g = zeros(nMom, nCells + 1);

if entropy_state.available
    left_alpha = entropy_state.left_alpha(:, 1);
    right_alpha = entropy_state.right_alpha(:, nCells);
else
    left_alpha = [];
    right_alpha = [];
end

[g(:, 1), ~] = boundary_interface_flux('left', uL(:, 1), model, quad_flux, opt_cfg, phys, t_stage, left_alpha, boundary_cache.left);

if use_par_interfaces && nCells > 1
    g_int = zeros(nMom, nCells - 1);
    parfor iif = 1:(nCells - 1)
        flux_cfg_int = struct('opt_cfg', opt_cfg);
        if entropy_state.available
            flux_cfg_int.alpha_left = entropy_state.right_alpha(:, iif);
            flux_cfg_int.alpha_right = entropy_state.left_alpha(:, iif + 1);
        end
        [g_int(:, iif), ~] = mm_kinetic_flux(uR(:, iif), uL(:, iif + 1), model, quad_flux, flux_cfg_int);
    end
    g(:, 2:nCells) = g_int;
else
    for iif = 1:(nCells - 1)
        flux_cfg_int = struct('opt_cfg', opt_cfg);
        if entropy_state.available
            flux_cfg_int.alpha_left = entropy_state.right_alpha(:, iif);
            flux_cfg_int.alpha_right = entropy_state.left_alpha(:, iif + 1);
        end
        [g(:, iif + 1), ~] = mm_kinetic_flux(uR(:, iif), uL(:, iif + 1), model, quad_flux, flux_cfg_int);
    end
end

[g(:, nCells + 1), ~] = boundary_interface_flux('right', uR(:, nCells), model, quad_flux, opt_cfg, phys, t_stage, right_alpha, boundary_cache.right);
end

function [stateL, stateR, rec_state] = prepare_entropy_reconstructed_states(u_stage, uL_in, uR_in, model, quad_flux, opt_cfg, stage_cache)
[nMom, nCells] = size(u_stage);
stateL = uL_in;
stateR = uR_in;

rec_state = struct();
rec_state.available = true;
rec_state.left_alpha = zeros(nMom, nCells);
rec_state.right_alpha = zeros(nMom, nCells);
rec_state.failed_cells = false(1, nCells);

for ic = 1:nCells
    warm_cache = get_reconstructed_warm_cache(stage_cache, ic);
    [alphaL, infoL] = mm_entropy_dual_solve(stateL(:, ic), model, quad_flux, opt_cfg, warm_cache);
    [alphaR, infoR] = mm_entropy_dual_solve(stateR(:, ic), model, quad_flux, opt_cfg, warm_cache);
    rec_state.left_alpha(:, ic) = alphaL(:);
    rec_state.right_alpha(:, ic) = alphaR(:);
    rec_state.failed_cells(ic) = ~(infoL.converged && infoR.converged);
end

if isfield(model, 'reflection_matrix') && ~isempty(model.reflection_matrix) && ...
        stage_is_mirrored(u_stage, model.reflection_matrix)
    for ic = 1:nCells
        jc = nCells + 1 - ic;
        if rec_state.failed_cells(ic) || rec_state.failed_cells(jc)
            rec_state.failed_cells(ic) = true;
            rec_state.failed_cells(jc) = true;
        end
    end
end

failed_idx = find(rec_state.failed_cells);
for k = 1:numel(failed_idx)
    ic = failed_idx(k);
    stateL(:, ic) = u_stage(:, ic);
    stateR(:, ic) = u_stage(:, ic);
    exact_cell_cache = get_exact_cell_cache(stage_cache, ic);
    alphaBar = get_field_or(exact_cell_cache, 'alpha', []);
    if isempty(alphaBar)
        [alphaBar, ~] = mm_entropy_dual_solve(u_stage(:, ic), model, quad_flux, opt_cfg, exact_cell_cache);
    end
    rec_state.left_alpha(:, ic) = alphaBar(:);
    rec_state.right_alpha(:, ic) = alphaBar(:);
end

rec_state.n_failed_cells = numel(failed_idx);
end

function [stage_cache, jac_state, timing] = prepare_stage_cell_cache(u_stage, model, quad_flux, opt_cfg, warm_stage_cache, use_char_recon, use_par_cells)
[nMom, nCells] = size(u_stage);
tPrep = tic;

cell_cache = cell(1, nCells);
info_cells = cell(1, nCells);
converged_cells = true(1, nCells);
iterations_cells = zeros(1, nCells);
flux_cell = zeros(nMom, nCells);
alpha_cells = zeros(nMom, nCells);
char_times = zeros(1, nCells);
entropy_times = zeros(1, nCells);

const_char = use_char_recon && ~model.needs_entropy && ...
    isfield(quad_flux, 'linear_char') && quad_flux.linear_char.success;

if const_char
    jac_state = build_constant_jacobian_state(quad_flux.linear_char);
else
    jac_state = struct();
    if use_char_recon
        jac_state.V = cell(1, nCells);
        jac_state.Vinv = cell(1, nCells);
        jac_state.lambda = cell(1, nCells);
    end
end

if model.needs_entropy
    if use_par_cells
        V_cells = cell(1, nCells);
        Vinv_cells = cell(1, nCells);
        lambda_cells = cell(1, nCells);
        parfor ic = 1:nCells
            warm_cell = get_warm_cell_cache(warm_stage_cache, ic);

            tEntropy = tic;
            [f_i, f_state] = mm_eval_flux_function(u_stage(:, ic), model, quad_flux, opt_cfg, warm_cell);
            entropy_times(ic) = toc(tEntropy);

            alpha_i = f_state.alpha(:);
            info_i = get_field_or(f_state.aux, 'info', struct('converged', true, 'iterations', 0));
            flux_cell(:, ic) = f_i(:);
            alpha_cells(:, ic) = alpha_i;
            info_cells{ic} = info_i;
            converged_cells(ic) = get_field_or(info_i, 'converged', true);
            iterations_cells(ic) = get_field_or(info_i, 'iterations', 0);

            cell_i = struct('alpha', alpha_i, 'flux', f_i(:), 'info', info_i);
            if use_char_recon
                tChar = tic;
                [V_i, Vinv_i, char_state_i] = mm_characteristic_basis(u_stage(:, ic), model, quad_flux, opt_cfg, cell_i);
                char_times(ic) = toc(tChar);
                V_cells{ic} = V_i;
                Vinv_cells{ic} = Vinv_i;
                lambda_cells{ic} = get_field_or(char_state_i, 'lambda', zeros(nMom, 1));
                cell_i.V = V_i;
                cell_i.Vinv = Vinv_i;
                cell_i.lambda = lambda_cells{ic};
            end
            cell_cache{ic} = cell_i;
        end

        if use_char_recon
            jac_state.V = V_cells;
            jac_state.Vinv = Vinv_cells;
            jac_state.lambda = lambda_cells;
        end
    else
        for ic = 1:nCells
            warm_cell = get_warm_cell_cache(warm_stage_cache, ic);

            tEntropy = tic;
            [f_i, f_state] = mm_eval_flux_function(u_stage(:, ic), model, quad_flux, opt_cfg, warm_cell);
            entropy_times(ic) = toc(tEntropy);

            alpha_i = f_state.alpha(:);
            info_i = get_field_or(f_state.aux, 'info', struct('converged', true, 'iterations', 0));
            flux_cell(:, ic) = f_i(:);
            alpha_cells(:, ic) = alpha_i;
            info_cells{ic} = info_i;
            converged_cells(ic) = get_field_or(info_i, 'converged', true);
            iterations_cells(ic) = get_field_or(info_i, 'iterations', 0);

            cell_i = struct('alpha', alpha_i, 'flux', f_i(:), 'info', info_i);
            if use_char_recon
                tChar = tic;
                [V_i, Vinv_i, char_state_i] = mm_characteristic_basis(u_stage(:, ic), model, quad_flux, opt_cfg, cell_i);
                char_times(ic) = toc(tChar);
                jac_state.V{ic} = V_i;
                jac_state.Vinv{ic} = Vinv_i;
                jac_state.lambda{ic} = get_field_or(char_state_i, 'lambda', zeros(nMom, 1));
                cell_i.V = V_i;
                cell_i.Vinv = Vinv_i;
                cell_i.lambda = jac_state.lambda{ic};
            end
            cell_cache{ic} = cell_i;
        end
    end
else
    alpha_cells = model.mass_matrix \ u_stage;
    if isfield(quad_flux, 'linear_char') && isfield(quad_flux.linear_char, 'K_flux') && ~isempty(quad_flux.linear_char.K_flux)
        flux_cell = quad_flux.linear_char.K_flux * alpha_cells;
    else
        psi_cell = quad_flux.B * alpha_cells;
        flux_cell = quad_flux.B.' * bsxfun(@times, quad_flux.w(:) .* quad_flux.mu(:), psi_cell);
    end
    for ic = 1:nCells
        info_i = struct('converged', true, 'iterations', 0, 'cached_alpha', alpha_cells(:, ic));
        info_cells{ic} = info_i;
        cell_cache{ic} = struct('alpha', alpha_cells(:, ic), 'flux', flux_cell(:, ic), 'info', info_i);
    end
end

stage_cache = struct();
stage_cache.alpha_cells = alpha_cells;
stage_cache.flux_cell = flux_cell;
stage_cache.cell_cache = cell_cache;
stage_cache.info_cells = info_cells;
stage_cache.converged_cells = converged_cells;
stage_cache.iterations_cells = iterations_cells;
stage_cache.linear_char = get_field_or(quad_flux, 'linear_char', struct());
stage_cache.jac_state = jac_state;

timing = struct();
timing.prepare_s = toc(tPrep);
timing.entropy_cell_solves_s = sum(entropy_times);
timing.characteristic_basis_s = sum(char_times);
end

function cache_out = extract_flux_cache_state(cache_in)
cache_out = struct();
if ~isstruct(cache_in) || isempty(fieldnames(cache_in))
    return;
end

if isfield(cache_in, 'cache_state') && isstruct(cache_in.cache_state)
    cache_out = cache_in.cache_state;
    return;
end

if isfield(cache_in, 'last_flux') && isstruct(cache_in.last_flux) && ...
        isfield(cache_in.last_flux, 'cache_state') && isstruct(cache_in.last_flux.cache_state)
    cache_out = cache_in.last_flux.cache_state;
    return;
end

cache_out = cache_in;
end

function warm_cache = get_stage_warm_cache(base_cache, prev_stage_cache)
warm_cache = struct();
if nargin >= 2 && isstruct(prev_stage_cache) && isfield(prev_stage_cache, 'alpha_cells') && ~isempty(prev_stage_cache.alpha_cells)
    warm_cache.alpha_cells = prev_stage_cache.alpha_cells;
    warm_cache.info_cells = get_field_or(prev_stage_cache, 'info_cells', cell(1, size(prev_stage_cache.alpha_cells, 2)));
else
    if isfield(base_cache, 'stage2_alpha_cells') && ~isempty(base_cache.stage2_alpha_cells)
        warm_cache.alpha_cells = base_cache.stage2_alpha_cells;
        warm_cache.info_cells = get_field_or(base_cache, 'stage2_info_cells', cell(1, size(base_cache.stage2_alpha_cells, 2)));
    elseif isfield(base_cache, 'last_alpha_cells') && ~isempty(base_cache.last_alpha_cells)
        warm_cache.alpha_cells = base_cache.last_alpha_cells;
        warm_cache.info_cells = get_field_or(base_cache, 'last_info_cells', cell(1, size(base_cache.last_alpha_cells, 2)));
    end
end

if isfield(base_cache, 'boundary_signature_cache')
    warm_cache.boundary_signature_cache = base_cache.boundary_signature_cache;
end
if isfield(base_cache, 'linear_char')
    warm_cache.linear_char = base_cache.linear_char;
end
end

function next_cache = build_next_flux_cache(base_cache, stage1_cache, stage2_cache)
next_cache = struct();
next_cache.last_alpha_cells = get_field_or(stage2_cache, 'alpha_cells', []);
next_cache.stage1_alpha_cells = get_field_or(stage1_cache, 'alpha_cells', []);
next_cache.stage2_alpha_cells = get_field_or(stage2_cache, 'alpha_cells', []);
next_cache.last_info_cells = get_field_or(stage2_cache, 'info_cells', {});
next_cache.stage1_info_cells = get_field_or(stage1_cache, 'info_cells', {});
next_cache.stage2_info_cells = get_field_or(stage2_cache, 'info_cells', {});
next_cache.last_converged_cells = get_field_or(stage2_cache, 'converged_cells', []);
next_cache.last_iterations_cells = get_field_or(stage2_cache, 'iterations_cells', []);
next_cache.linear_char = get_field_or(stage2_cache, 'linear_char', get_field_or(base_cache, 'linear_char', struct()));
next_cache.boundary_signature_cache = get_field_or(base_cache, 'boundary_signature_cache', struct());
end

function warm_cell = get_warm_cell_cache(stage_cache, ic)
warm_cell = struct();
if ~isstruct(stage_cache)
    return;
end

alpha0 = [];
if isfield(stage_cache, 'alpha_cells') && size(stage_cache.alpha_cells, 2) >= ic
    alpha0 = stage_cache.alpha_cells(:, ic);
elseif isfield(stage_cache, 'last_alpha_cells') && size(stage_cache.last_alpha_cells, 2) >= ic
    alpha0 = stage_cache.last_alpha_cells(:, ic);
end

if ~isempty(alpha0)
    warm_cell.last_alpha = alpha0(:);
    warm_cell.alpha0 = alpha0(:);
end

if isfield(stage_cache, 'info_cells') && numel(stage_cache.info_cells) >= ic && ~isempty(stage_cache.info_cells{ic})
    warm_cell.info = stage_cache.info_cells{ic};
end
end

function warm_cell = get_reconstructed_warm_cache(stage_cache, ic)
warm_cell = get_warm_cell_cache(stage_cache, ic);
exact_cell = get_exact_cell_cache(stage_cache, ic);
if isfield(exact_cell, 'alpha') && ~isempty(exact_cell.alpha)
    warm_cell.last_alpha = exact_cell.alpha(:);
    warm_cell.alpha0 = exact_cell.alpha(:);
end
end

function exact_cell = get_exact_cell_cache(stage_cache, ic)
exact_cell = struct();
if ~isstruct(stage_cache)
    return;
end

if isfield(stage_cache, 'cell_cache') && numel(stage_cache.cell_cache) >= ic && ~isempty(stage_cache.cell_cache{ic})
    exact_cell = stage_cache.cell_cache{ic};
end

if ~isfield(exact_cell, 'alpha') && isfield(stage_cache, 'alpha_cells') && size(stage_cache.alpha_cells, 2) >= ic
    exact_cell.alpha = stage_cache.alpha_cells(:, ic);
end
if ~isfield(exact_cell, 'flux') && isfield(stage_cache, 'flux_cell') && size(stage_cache.flux_cell, 2) >= ic
    exact_cell.flux = stage_cache.flux_cell(:, ic);
end
if ~isfield(exact_cell, 'info') && isfield(stage_cache, 'info_cells') && numel(stage_cache.info_cells) >= ic
    exact_cell.info = stage_cache.info_cells{ic};
end
end

function jac_state = build_constant_jacobian_state(linear_char)
jac_state = struct();
jac_state.constant = true;
jac_state.V_const = linear_char.V;
jac_state.Vinv_const = linear_char.Vinv;
jac_state.lambda_const = linear_char.lambda;
jac_state.method = get_field_or(linear_char, 'method', 'linear-generalized');
end

function [flux, uB] = boundary_interface_flux(side, uCell, model, quad_flux, opt_cfg, phys, t, cell_alpha, side_cache)
if nargin < 8
    cell_alpha = [];
end
if nargin < 9 || isempty(side_cache)
    side_cache = prepare_boundary_side_cache(side, phys, t, model, quad_flux);
end

uB = side_cache.uB;

if ~side_cache.has_psi
    flux_cfg = struct('opt_cfg', opt_cfg);
    if isfield(side_cache, 'alpha') && ~isempty(side_cache.alpha)
        if strcmp(side, 'left')
            flux_cfg.alpha_left = side_cache.alpha(:);
        else
            flux_cfg.alpha_right = side_cache.alpha(:);
        end
    end
    if ~isempty(cell_alpha)
        if strcmp(side, 'left')
            flux_cfg.alpha_right = cell_alpha(:);
        else
            flux_cfg.alpha_left = cell_alpha(:);
        end
    end

    if strcmp(side, 'left')
        [flux, ~] = mm_kinetic_flux(uB, uCell, model, quad_flux, flux_cfg);
    else
        [flux, ~] = mm_kinetic_flux(uCell, uB, model, quad_flux, flux_cfg);
    end
    return;
end

if ~isempty(cell_alpha)
    if model.needs_entropy
        psiCell = exp(min(quad_flux.B * cell_alpha(:), 700));
    else
        psiCell = quad_flux.B * cell_alpha(:);
    end
else
    [psiCell, ~, ~] = mm_eval_ansatz(uCell, model, quad_flux, opt_cfg, struct());
end
psiCell = max(psiCell(:), 0.0);

if strcmp(side, 'left')
    psiOutMinus = psiCell(quad_flux.mu <= 0);
    flux = quad_flux.B_plus.' * (quad_flux.w_plus .* quad_flux.mu_plus .* side_cache.psi_plus) + ...
           quad_flux.B_minus.' * (quad_flux.w_minus .* quad_flux.mu_minus .* psiOutMinus);
else
    psiOutPlus = psiCell(quad_flux.mu >= 0);
    flux = quad_flux.B_plus.' * (quad_flux.w_plus .* quad_flux.mu_plus .* psiOutPlus) + ...
           quad_flux.B_minus.' * (quad_flux.w_minus .* quad_flux.mu_minus .* side_cache.psi_minus);
end
end

function boundary_cache = prepare_boundary_cache(phys, t, model, quad_flux)
boundary_cache = struct();
boundary_cache.left = prepare_boundary_side_cache('left', phys, t, model, quad_flux);
boundary_cache.right = prepare_boundary_side_cache('right', phys, t, model, quad_flux);
end

function side_cache = prepare_boundary_side_cache(side, phys, t, model, quad_flux)
rhoB = boundary_density(side, t, phys);
side_cache = struct();
side_cache.rho = rhoB;
side_cache.has_psi = has_boundary_psi(side, phys);

if side_cache.has_psi
    psi_full = eval_boundary_psi(side, quad_flux.mu, t, phys, rhoB);
    psi_full = max(psi_full(:), 0.0);
    side_cache.psi_full = psi_full;
    side_cache.uB = quad_flux.B.' * (quad_flux.w .* psi_full);
    side_cache.fB = quad_flux.B.' * (quad_flux.w .* quad_flux.mu .* psi_full);
    side_cache.psi_plus = psi_full(quad_flux.mu >= 0);
    side_cache.psi_minus = psi_full(quad_flux.mu <= 0);
else
    psi_iso = (rhoB / model.h1) * ones(size(quad_flux.mu));
    side_cache.psi_full = psi_iso(:);
    side_cache.psi_plus = psi_iso(quad_flux.mu >= 0);
    side_cache.psi_minus = psi_iso(quad_flux.mu <= 0);
    side_cache.uB = model.b_iso * rhoB;
    side_cache.fB = quad_flux.B.' * (quad_flux.w .* quad_flux.mu .* psi_iso(:));
    if model.needs_entropy
        side_cache.alpha = model.alpha1 * log(max(rhoB, realmin) / model.h1);
    else
        side_cache.alpha = model.mass_matrix \ side_cache.uB;
    end
end
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
timing.stage_cache_preparation_s = get_timing_field(s1, 'stage_cache_preparation_s') + get_timing_field(s2, 'stage_cache_preparation_s');
timing.entropy_cell_solves_s = get_timing_field(s1, 'entropy_cell_solves_s') + get_timing_field(s2, 'entropy_cell_solves_s');
timing.jacobian_s = get_timing_field(s1, 'jacobian_s') + get_timing_field(s2, 'jacobian_s');
timing.characteristic_basis_s = get_timing_field(s1, 'characteristic_basis_s') + get_timing_field(s2, 'characteristic_basis_s');
timing.reconstruct_s = get_timing_field(s1, 'reconstruct_s') + get_timing_field(s2, 'reconstruct_s');
timing.entropy_reconstructed_solves_s = get_timing_field(s1, 'entropy_reconstructed_solves_s') + get_timing_field(s2, 'entropy_reconstructed_solves_s');
timing.high_order_flux_s = get_timing_field(s1, 'high_order_flux_s') + get_timing_field(s2, 'high_order_flux_s');
timing.low_order_flux_s = get_timing_field(s1, 'low_order_flux_s') + get_timing_field(s2, 'low_order_flux_s');
timing.mcl_bisection_s = get_timing_field(s1, 'mcl_bisection_s') + get_timing_field(s2, 'mcl_bisection_s');
timing.limiter_and_flux_s = get_timing_field(s1, 'limiter_and_flux_s') + get_timing_field(s2, 'limiter_and_flux_s');
timing.rhs_build_s = get_timing_field(s1, 'rhs_build_s') + get_timing_field(s2, 'rhs_build_s');
end

function v = get_timing_field(stageState, name)
v = 0.0;
if isfield(stageState, 'timing') && isfield(stageState.timing, name)
    v = stageState.timing.(name);
end
end
