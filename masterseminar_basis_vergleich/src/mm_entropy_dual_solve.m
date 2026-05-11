function [alpha, info] = mm_entropy_dual_solve(u, model, quad, opt_cfg, cache_state)
%MM_ENTROPY_DUAL_SOLVE Newton/Armijo solver for entropy dual problem.
% References (Seminarquelle 2):
% - Dual objective, gradient, Hessian: Eq. (4.19)-(4.21)
% - Density normalization / multiplier rescaling: Eq. (4.22)-(4.23)
% - Stopping criteria and practical safeguards: Eq. (4.25), (5.4), (5.5)

if nargin < 5 || isempty(cache_state)
    cache_state = struct();
end

u = u(:);
n = model.nMom;

if ~model.needs_entropy
    alpha = model.mass_matrix \ u;
    info = struct('converged', true, 'iterations', 0, 'regularization_r', 0, ...
        'used_linear_closure', true, 'vacuum_regularized', false, ...
        'cached_alpha', alpha, 'used_isotropic_fallback', false);
    return;
end

opt = fill_defaults(opt_cfg, model);

if strcmp(model.family, 'partial')
    [alpha, info] = solve_partial_entropy_blocks(u, model, quad, opt, cache_state);
    return;
end

check_fun = @(v) mm_is_realizable(v, model, quad, struct('epsR', 0));
[alpha, info] = solve_entropy_problem(u, model, quad.B, quad.w, opt, cache_state, check_fun);

end

function [alpha, info] = solve_partial_entropy_blocks(u, model, quad, opt, cache_state)
k = model.kIntervals;
alpha = zeros(model.nMom, 1);
info = init_info();
info.method = 'partial-block';
info.block_info = cell(1, k);

alpha0 = get_cached_alpha0(cache_state, model.nMom);

all_converged = true;
max_iter = 0;
max_reg = 0.0;
any_vac = false;
all_crit1 = true;
all_crit2 = true;
any_iso_fallback = false;

for ib = 1:k
    idx = (2 * ib - 1):(2 * ib);
    block_model = build_partial_block_model(model, ib);
    block_quad.B = quad.local_B_by_interval{ib};
    block_quad.w = quad.local_w_by_interval{ib};

    block_cache = struct();
    if numel(alpha0) == model.nMom
        block_cache.last_alpha = alpha0(idx);
        block_cache.alpha0 = alpha0(idx);
    end

    a = block_model.mu_edges(1);
    b = block_model.mu_edges(2);
    check_fun = @(v) check_partial_block_realizable(v, a, b);
    [alpha_b, info_b] = solve_entropy_problem(u(idx), block_model, block_quad.B, block_quad.w, opt, block_cache, check_fun);

    alpha(idx) = alpha_b(:);
    info.block_info{ib} = info_b;
    all_converged = all_converged && info_b.converged;
    max_iter = max(max_iter, info_b.iterations);
    max_reg = max(max_reg, info_b.regularization_r);
    any_vac = any_vac || info_b.vacuum_regularized;
    all_crit1 = all_crit1 && info_b.criterion1;
    all_crit2 = all_crit2 && info_b.criterion2;
    any_iso_fallback = any_iso_fallback || info_b.used_isotropic_fallback;
end

info.converged = all_converged;
info.iterations = max_iter;
info.regularization_r = max_reg;
info.vacuum_regularized = any_vac;
info.criterion1 = all_crit1;
info.criterion2 = all_crit2;
info.used_isotropic_fallback = any_iso_fallback;
info.cached_alpha = alpha;
end

function [alpha, info] = solve_entropy_problem(u, model, B, w, opt, cache_state, check_fun)
u = u(:);
B = B;
w = w(:);
n = model.nMom;

info = init_info();

rho = model.alpha1.' * u;
if rho <= 0
    rho = opt.rho_vac;
    u = model.u_iso_vec * rho;
    info.vacuum_regularized = true;
end

if rho < opt.rho_vac
    rho = opt.rho_vac;
    u = model.u_iso_vec * rho;
    info.vacuum_regularized = true;
end

beta0 = model.alpha1 * log(1.0 / model.h1);
alpha_cache = get_cached_alpha0(cache_state, n);
if numel(alpha_cache) == n
    beta0 = alpha_cache(:);
end

for rr = 1:numel(opt.regularization_r)
    r = opt.regularization_r(rr);
    u_try = (1 - r) * u + r * (model.G * u);
    rho_try = model.alpha1.' * u_try;
    phi = u_try / rho_try;

    beta = beta0;
    d_last = zeros(n, 1);

    for k = 1:opt.max_iter
        z = B * beta;
        ez = safe_exp(z);

        u_beta = B.' * (w .* ez);
        g = u_beta - phi;
        gnorm = norm(g, 2);

        tau0 = compute_tau0(phi, rho_try, opt.tau, model);
        crit1 = gnorm < min(tau0, opt.tau);

        if crit1
            varrho = sum(w .* ez);
            alpha_try = beta + model.alpha1 * log(rho_try / varrho);
            u_alpha = B.' * (w .* safe_exp(B * alpha_try));

            need_check = true;
            if ~model.is_hat
                need_check = (1 - opt.eps_gamma) < exp(-(norm(d_last, 1) + abs(log(max(varrho, eps)))));
            end

            crit2 = true;
            if need_check
                crit2 = check_fun(u_try - (1 - opt.eps_gamma) * u_alpha);
            end

            info.criterion1 = crit1;
            info.criterion2 = crit2;

            if crit2
                alpha = alpha_try;
                info.converged = true;
                info.iterations = k;
                info.regularization_r = r;
                info.cached_alpha = alpha;
                return;
            end
        end

        if gnorm < opt.grad_tol
            break;
        end

        H = weighted_basis_gram(B, w .* ez);
        d = compute_newton_direction(H, g, opt.use_change_of_basis);
        if any(~isfinite(d))
            break;
        end

        p0 = sum(w .* ez) - phi.' * beta;
        step = 1.0;
        accepted = false;

        while step >= opt.step_min
            bt = beta + step * d;
            zt = B * bt;
            ezt = safe_exp(zt);
            p1 = sum(w .* ezt) - phi.' * bt;

            if p1 < p0 + opt.armijo_c * step * (g.' * d)
                accepted = true;
                beta = bt;
                d_last = d;
                break;
            end
            step = step * opt.armijo_backtrack;
        end

        if ~accepted
            break;
        end

        info.iterations = max(info.iterations, k);
    end
end

rho_fallback = max(model.alpha1.' * u, opt.rho_vac);
alpha = model.alpha1 * log(rho_fallback / model.h1);
info.converged = false;
info.used_isotropic_fallback = true;
info.cached_alpha = alpha;
end

function info = init_info()
info = struct();
info.converged = false;
info.iterations = 0;
info.regularization_r = NaN;
info.vacuum_regularized = false;
info.criterion1 = false;
info.criterion2 = false;
info.used_isotropic_fallback = false;
info.used_linear_closure = false;
info.cached_alpha = [];
end

function d = compute_newton_direction(H, g, use_change_of_basis)
reg = 1.0e-12;
n = size(H, 1);

if use_change_of_basis
    dH = full(diag(H));
    s = 1 ./ sqrt(max(dH, 1.0e-30));
    Hs = H .* (s * s.');
    gs = s .* g;
    Hs = Hs + reg * speye(n);
    ds = solve_linear(Hs, -gs, reg);
    d = s .* ds;
else
    H = H + reg * speye(n);
    d = solve_linear(H, -g, reg);
end
end

function x = solve_linear(A, b, reg)
n = size(A, 1);
if size(A, 2) ~= n
    x = A \ b;
    if any(~isfinite(x))
        x = pinv(full(A)) * b;
    end
    return;
end

[L, p] = chol(A, 'lower');
if p == 0
    x = L' \ (L \ b);
else
    regA = max(reg, 1.0e-10 * max(1.0, norm(A, 1)));
    x = (A + regA * speye(n)) \ b;
end

if any(~isfinite(x))
    x = pinv(full(A)) * b;
end
end

function G = weighted_basis_gram(B, weights)
weights = weights(:);
if issparse(B)
    G = B.' * spdiags(weights, 0, numel(weights), numel(weights)) * B;
else
    G = B.' * bsxfun(@times, weights, B);
end
G = 0.5 * (G + G.');
end

function tf = check_partial_block_realizable(u, a, b)
u = u(:);
u0 = u(1);
u1 = u(2);
tf = (u0 > 0) && ((u1 - a * u0) > 0) && ((b * u0 - u1) > 0);
end

function block_model = build_partial_block_model(model, ib)
a = model.mu_edges(ib);
b = model.mu_edges(ib + 1);
h = b - a;
m1 = 0.5 * (b^2 - a^2) / h;

block_model = struct();
block_model.nMom = 2;
block_model.alpha1 = [1.0; 0.0];
block_model.h1 = h;
block_model.u_iso_vec = [1.0; m1];
block_model.G = block_model.u_iso_vec * block_model.alpha1.';
block_model.is_hat = false;
block_model.family = 'partial';
block_model.mu_edges = [a; b];
end

function alpha0 = get_cached_alpha0(cache_state, n)
alpha0 = [];
if isfield(cache_state, 'last_alpha') && numel(cache_state.last_alpha) == n
    alpha0 = cache_state.last_alpha(:);
elseif isfield(cache_state, 'alpha0') && numel(cache_state.alpha0) == n
    alpha0 = cache_state.alpha0(:);
elseif isfield(cache_state, 'alpha') && numel(cache_state.alpha) == n
    alpha0 = cache_state.alpha(:);
elseif isfield(cache_state, 'alpha_cells') && size(cache_state.alpha_cells, 2) == 1 && size(cache_state.alpha_cells, 1) == n
    alpha0 = cache_state.alpha_cells(:, 1);
elseif isfield(cache_state, 'last_alpha_cells') && size(cache_state.last_alpha_cells, 2) == 1 && size(cache_state.last_alpha_cells, 1) == n
    alpha0 = cache_state.last_alpha_cells(:, 1);
end
end

function tau0 = compute_tau0(phi, rho, tau, model)
n = model.nMom;
if strcmp(model.family, 'full')
    tau0 = tau / ((1 + norm(phi, 2)^2) * rho + tau);
else
    tau0 = tau / ((1 + sqrt(n) * norm(phi, 2)^2) * rho + sqrt(n) * tau);
end
end

function y = safe_exp(x)
y = exp(min(x, 700));
end

function opt = fill_defaults(opt_cfg, model)
opt = struct();
opt.tau = get_field_or(opt_cfg, 'tau', 1.0e-9);
opt.eps_gamma = get_field_or(opt_cfg, 'eps_gamma', 1.0e-2);
opt.max_iter = get_field_or(opt_cfg, 'max_iter', 200);
opt.armijo_c = get_field_or(opt_cfg, 'armijo_c', 1.0e-3);
opt.armijo_backtrack = get_field_or(opt_cfg, 'armijo_backtrack', 0.5);
opt.step_min = get_field_or(opt_cfg, 'step_min', 1.0e-12);
opt.grad_tol = get_field_or(opt_cfg, 'grad_tol', 1.0e-11);
opt.rho_vac = get_field_or(opt_cfg, 'rho_vac', 1.0e-8);
opt.regularization_r = get_field_or(opt_cfg, 'regularization_r', [0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1]);

if isfield(opt_cfg, 'use_change_of_basis')
    opt.use_change_of_basis = logical(opt_cfg.use_change_of_basis);
else
    opt.use_change_of_basis = model.use_change_of_basis_default;
end

if model.is_partial && model.needs_entropy && isfield(opt_cfg, 'use_change_of_basis_partial_entropy')
    opt.use_change_of_basis = logical(opt_cfg.use_change_of_basis_partial_entropy);
end
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
