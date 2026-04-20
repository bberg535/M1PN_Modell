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
        'used_linear_closure', true, 'vacuum_regularized', false);
    return;
end

opt = fill_defaults(opt_cfg, model);

info = struct();
info.converged = false;
info.iterations = 0;
info.regularization_r = NaN;
info.vacuum_regularized = false;
info.criterion1 = false;
info.criterion2 = false;

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

beta0 = zeros(n, 1);
if isfield(cache_state, 'last_alpha') && numel(cache_state.last_alpha) == n
    beta0 = cache_state.last_alpha(:);
elseif isfield(cache_state, 'alpha0') && numel(cache_state.alpha0) == n
    beta0 = cache_state.alpha0(:);
end

B = quad.B;
w = quad.w;

for rr = 1:numel(opt.regularization_r)
    r = opt.regularization_r(rr);
    u_try = (1 - r) * u + r * (model.G * u);
    rho_try = model.alpha1.' * u_try;
    phi = u_try / rho_try;

    beta = beta0;
    d_last = zeros(n, 1);

    for k = 1:opt.max_iter
        % u(beta)=<b exp(b·beta)> and g(beta)=u(beta)-phi, cf. Eq. (4.20).
        z = B * beta;
        ez = safe_exp(z);

        u_beta = B.' * (w .* ez);
        g = u_beta - phi;
        gnorm = norm(g, 2);

        tau0 = compute_tau0(phi, rho_try, opt.tau, model);
        % First stopping criterion (Eq. (4.25a) with practical tau0 from Section 5.1).
        crit1 = gnorm < min(tau0, opt.tau);

        if crit1
            varrho = sum(w .* ez);
            % Rescaling of multipliers to preserve local density, Eq. (4.23).
            alpha_try = beta + model.alpha1 * log(rho_try / varrho);
            u_alpha = B.' * (w .* safe_exp(B * alpha_try));

            need_check = true;
            if ~model.is_hat
                need_check = (1 - opt.eps_gamma) < exp(-(norm(d_last, 1) + abs(log(max(varrho, eps)))));
            end

            crit2 = true;
            if need_check
                % Practical realizability safeguard, Eq. (5.5) variant.
                [crit2, ~] = mm_is_realizable(u_try - (1 - opt.eps_gamma) * u_alpha, model, quad, struct('epsR', 0));
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

        % Hessian H=<b b^T exp(b·beta)>, Eq. (4.21).
        H = B.' * (((w .* ez) .* B));

        d = compute_newton_direction(H, g, opt.use_change_of_basis);
        if any(~isfinite(d))
            break;
        end

        % Dual objective p(beta), Eq. (4.19).
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

% Fallback to isotropic multiplier.
rho_fallback = max(model.alpha1.' * u, opt.rho_vac);
alpha = model.alpha1 * log(rho_fallback / model.h1);
info.converged = false;
info.cached_alpha = alpha;

end

function d = compute_newton_direction(H, g, use_change_of_basis)

reg = 1.0e-12;
n = size(H, 1);

if use_change_of_basis
    dH = diag(H);
    s = 1 ./ sqrt(max(dH, 1e-30));
    S = diag(s);
    Hs = S * H * S;
    gs = S * g;
    Hs = Hs + reg * eye(n);
    ds = solve_linear(Hs, -gs);
    d = S * ds;
else
    H = H + reg * eye(n);
    d = solve_linear(H, -g);
end
end

function x = solve_linear(A, b)
n = size(A, 1);
if size(A, 2) ~= n
    x = pinv(A) * b;
    return;
end

rA = rcond(A);
if ~isfinite(rA) || rA < 1.0e-12
    % Avoid warning storms from near-singular solves in long PMMn runs.
    regA = max(1.0e-12, 1.0e-10 * norm(A, 1));
    Areg = A + regA * eye(n);
    rReg = rcond(Areg);
    if isfinite(rReg) && rReg >= 1.0e-12
        x = Areg \ b;
    else
        x = pinv(A) * b;
    end
else
    x = A \ b;
end

if any(~isfinite(x))
    x = pinv(A) * b;
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
