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
lim_type = normalize_limiter_type(get_field_or(lim_cfg, 'type', 'paper'));

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

        jac_state = struct();
        jac_state.V = cell(1, nCells);
        jac_state.Vinv = cell(1, nCells);
        tJac = tic;
        for ic = 1:nCells
            if rec_cfg.use_characteristic
                J = flux_jacobian_fd(u_stage(:, ic), model, quad_flux, opt_cfg);
                [V, D] = eig(J);
                if any(~isfinite(V(:))) || any(~isfinite(diag(D)))
                    V = eye(nMom);
                end
                jac_state.V{ic} = V;
                jac_state.Vinv{ic} = pinv(V);
            else
                jac_state.V{ic} = eye(nMom);
                jac_state.Vinv{ic} = eye(nMom);
            end
        end
        tJacobian = toc(tJac);

        tRec = tic;
        [uL_rec, uR_rec, rec_state] = mm_reconstruct_characteristic(u_stage, jac_state, grid, rec_cfg);
        tReconstruct = toc(tRec);

        tLimFlux = tic;
        if strcmp(lim_type, 'mcl')
            [g, lim_state] = mcl_limited_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg);
        else
            lim_cfg_apply = lim_cfg;
            if get_field_or(lim_cfg, 'paper_lp_characteristic', false)
                lim_cfg_apply.characteristic_basis = jac_state;
            end
            [uL_lim, uR_lim, lim_state] = mm_apply_realizability_limiter(u_stage, uL_rec, uR_rec, model, lim_cfg_apply, quad_flux);
            g = zeros(nMom, nCells + 1);

            % Left boundary interface.
            uBLeft = model.b_iso * boundary_density('left', t_stage, phys);
            [g(:, 1), ~] = mm_kinetic_flux(uBLeft, uL_lim(:, 1), model, quad_flux, struct('opt_cfg', opt_cfg));

            % Interior interfaces.
            for iif = 1:(nCells - 1)
                [g(:, iif + 1), ~] = mm_kinetic_flux(uR_lim(:, iif), uL_lim(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
            end

            % Right boundary interface.
            uBRight = model.b_iso * boundary_density('right', t_stage, phys);
            [g(:, nCells + 1), ~] = mm_kinetic_flux(uR_lim(:, nCells), uBRight, model, quad_flux, struct('opt_cfg', opt_cfg));
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

function [g, lim_state] = mcl_limited_flux(u_stage, uL_rec, uR_rec, model, phys, t_stage, quad_flux, opt_cfg, lim_cfg)
[nMom, nCells] = size(u_stage);
g_ho = zeros(nMom, nCells + 1);

uBLeft = model.b_iso * boundary_density('left', t_stage, phys);
uBRight = model.b_iso * boundary_density('right', t_stage, phys);

% High-order fluxes from reconstructed states.
[g_ho(:, 1), ~] = mm_kinetic_flux(uBLeft, uL_rec(:, 1), model, quad_flux, struct('opt_cfg', opt_cfg));
for iif = 1:(nCells - 1)
    [g_ho(:, iif + 1), ~] = mm_kinetic_flux(uR_rec(:, iif), uL_rec(:, iif + 1), model, quad_flux, struct('opt_cfg', opt_cfg));
end
[g_ho(:, nCells + 1), ~] = mm_kinetic_flux(uR_rec(:, nCells), uBRight, model, quad_flux, struct('opt_cfg', opt_cfg));

g_lo = zeros(nMom, nCells + 1);
g = zeros(nMom, nCells + 1);
theta = ones(1, nCells + 1);
bar_ok = true(1, nCells + 1);

lambdaMCL = get_field_or(lim_cfg, 'mcl_lambda', 1.0);
lambdaMCL = max(lambdaMCL, 1.0e-12);
maxBisect = get_field_or(lim_cfg, 'mcl_bisect_iter', 40);
epsR = get_field_or(lim_cfg, 'epsR', 0.0);

for iif = 1:(nCells + 1)
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

    g_lo(:, iif) = 0.5 * (fL + fR) - 0.5 * lambdaMCL * (UR - UL);
    Ubar = 0.5 * (UL + UR) - (fR - fL) / (2 * lambdaMCL);

    A = g_ho(:, iif) - g_lo(:, iif);

    [okBar, ~] = mm_is_realizable(Ubar, model, quad_flux, struct('epsR', epsR));
    if ~okBar
        theta(iif) = 0.0;
        bar_ok(iif) = false;
        g(:, iif) = g_lo(:, iif);
        continue;
    end

    % MCL-style convex limiting on antidiffusive interface flux.
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

    theta(iif) = lo;
    g(:, iif) = g_lo(:, iif) + theta(iif) * A;
end

lim_state = struct();
lim_state.theta = theta;
lim_state.mode = 'mcl';
lim_state.bar_state_feasible = bar_ok;
lim_state.g_low = g_lo;
lim_state.g_high = g_ho;

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
