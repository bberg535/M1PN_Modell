function [V, Vinv, state] = mm_characteristic_basis(u, model, quad_flux, opt_cfg, cache_state)
%MM_CHARACTERISTIC_BASIS Characteristic basis for flux reconstruction.
% Reference: Seminarquelle 2, Section 5.2.
% Computes the right eigenvectors of f'(u)=J(u)H(u)^{-1} through the
% symmetric generalized eigenproblem J ztilde = lambda H ztilde.

if nargin < 5 || isempty(cache_state)
    cache_state = struct();
end

n = model.nMom;
need_inverse = (nargout >= 2);
opt = characteristic_options(opt_cfg);

V = eye(n);
Vinv = eye(n);
state = struct('success', false, 'method', 'identity', 'lambda', zeros(n, 1), ...
    'residual_inf', NaN, 'condition_est', NaN, 'constant', false);

[pre_basis, has_pre_basis] = extract_precomputed_basis(cache_state, quad_flux, model);
if has_pre_basis
    V = pre_basis.V;
    if need_inverse
        Vinv = pre_basis.Vinv;
    end
    state.success = true;
    state.method = get_field_or(pre_basis, 'method', 'precomputed');
    state.lambda = get_field_or(pre_basis, 'lambda', zeros(n, 1));
    state.constant = get_field_or(pre_basis, 'constant', false);

    if opt.check_basis_condition
        state.condition_est = rcond(V);
    end
    if opt.check_eigen_residual && isfield(pre_basis, 'J') && isfield(pre_basis, 'H') && isfield(pre_basis, 'Z')
        state.residual_inf = generalized_residual_norm(pre_basis.J, pre_basis.H, pre_basis.Z, state.lambda);
    end
    return;
end

try
    [J, H, matrix_state] = flux_jacobian_factors(u, model, quad_flux, opt_cfg, cache_state);
    state.matrix = matrix_state;

    if strcmp(get_field_or(matrix_state, 'kind', ''), 'partial-blocks')
        [Vtry, lambda, ok, Ztry, block_state] = partial_block_basis(matrix_state.J_blocks, matrix_state.H_blocks, opt);
        method = get_field_or(block_state, 'method', 'partial-block-explicit');
        state.block = block_state;
    else
        [Vtry, lambda, ok, Ztry] = generalized_flux_basis(J, H);
        method = 'symmetric-generalized';
    end

    if ~ok && opt.safe_characteristic_fallback && ~strcmp(get_field_or(matrix_state, 'kind', ''), 'partial-blocks')
        F = full(J / H);
        [Vtry, lambda, ok] = nonsymmetric_real_basis(F);
        Ztry = [];
        method = 'nonsymmetric-real-fallback';
    end

    if ~ok
        return;
    end

    [lambda, idx] = sort(real(lambda(:)), 'ascend');
    Vtry = real(Vtry(:, idx));
    if ~isempty(Ztry)
        Ztry = real(Ztry(:, idx));
    end

    [Vtry, Ztry, ok] = normalize_basis_columns(Vtry, Ztry);
    if ~ok
        return;
    end

    if opt.check_basis_condition
        state.condition_est = rcond(Vtry);
        if ~isfinite(state.condition_est) || state.condition_est < 1.0e-12
            return;
        end
    end

    if need_inverse
        VinvTry = Vtry \ eye(n);
        if any(~isfinite(VinvTry(:)))
            return;
        end
    else
        VinvTry = eye(n);
    end

    if opt.check_eigen_residual && ~isempty(Ztry)
        if strcmp(get_field_or(matrix_state, 'kind', ''), 'partial-blocks')
            residual = generalized_residual_block_norm(matrix_state.J_blocks, matrix_state.H_blocks, Ztry, lambda);
        else
            residual = generalized_residual_norm(J, H, Ztry, lambda);
        end
        state.residual_inf = residual;
        if ~isfinite(residual) || residual > opt.eigen_residual_tol
            return;
        end
    end

    V = Vtry;
    Vinv = VinvTry;
    state.success = true;
    state.method = method;
    state.lambda = lambda;
catch err
    state.error = err.message;
end

end

function [pre_basis, tf] = extract_precomputed_basis(cache_state, quad_flux, model)
pre_basis = struct();
tf = false;

if isstruct(cache_state)
    if get_field_or(cache_state, 'constant', false) && ...
            isfield(cache_state, 'V_const') && isfield(cache_state, 'Vinv_const')
        pre_basis = struct('V', cache_state.V_const, 'Vinv', cache_state.Vinv_const, ...
            'lambda', get_field_or(cache_state, 'lambda_const', zeros(model.nMom, 1)), ...
            'method', get_field_or(cache_state, 'method', 'cached-constant'), ...
            'constant', true);
        tf = true;
        return;
    end

    if isfield(cache_state, 'linear_char')
        pre_basis = cache_state.linear_char;
        tf = isfield(pre_basis, 'success') && pre_basis.success;
        if tf
            return;
        end
    end
end

if ~model.needs_entropy && isfield(quad_flux, 'linear_char')
    pre_basis = quad_flux.linear_char;
    tf = isfield(pre_basis, 'success') && pre_basis.success;
end
end

function [J, H, state] = flux_jacobian_factors(u, model, quad_flux, opt_cfg, cache_state)
B = quad_flux.B;
w = quad_flux.w(:);
mu = quad_flux.mu(:);

state = struct();
if model.needs_entropy
    [psi, alpha, aux] = mm_eval_ansatz(u, model, quad_flux, opt_cfg, cache_state);
    psi = max(real(psi(:)), realmin);
    weights = w .* psi;
    state.alpha = alpha;
    state.ansatz = aux;
    state.closure = 'entropy';
    state.psi = psi;
else
    weights = w;
    if isfield(cache_state, 'alpha') && numel(cache_state.alpha) == model.nMom
        state.alpha = cache_state.alpha(:);
    else
        state.alpha = model.mass_matrix \ u(:);
    end
    state.ansatz = struct('info', struct('converged', true, 'iterations', 0));
    state.closure = 'linear';
end

if strcmp(model.family, 'partial') && model.needs_entropy
    [Jblk, Hblk, block_moments] = partial_block_factors(weights, mu, quad_flux, model);
    J = sparse([]);
    H = sparse([]);
    state.kind = 'partial-blocks';
    state.J_blocks = Jblk;
    state.H_blocks = Hblk;
    state.block_moments = block_moments;
else
    J = weighted_basis_gram(B, weights .* mu);
    if model.needs_entropy
        H = weighted_basis_gram(B, weights);
    else
        H = model.mass_matrix;
    end
    state.kind = 'dense';
end
end

function [Jblk, Hblk, moments] = partial_block_factors(weights, mu, quad_flux, model)
k = model.kIntervals;
Jblk = zeros(2, 2, k);
Hblk = zeros(2, 2, k);
moments = zeros(4, k);

for ib = 1:k
    if isfield(quad_flux, 'interval_qidx') && numel(quad_flux.interval_qidx) >= ib
        qidx = quad_flux.interval_qidx{ib};
    else
        qidx = find(quad_flux.interval_id == ib);
    end

    wj = weights(qidx);
    muj = mu(qidx);
    m0 = sum(wj);
    m1 = sum(wj .* muj);
    m2 = sum(wj .* (muj.^2));
    m3 = sum(wj .* (muj.^3));

    Hblk(:, :, ib) = [m0, m1; m1, m2];
    Jblk(:, :, ib) = [m1, m2; m2, m3];
    moments(:, ib) = [m0; m1; m2; m3];
end
end

function [V, lambda, ok, Z, state] = partial_block_basis(J_blocks, H_blocks, opt)
k = size(J_blocks, 3);
n = 2 * k;
V = zeros(n, n);
Z = zeros(n, n);
lambda = zeros(n, 1);
ok = true;
methods = cell(1, k);

for ib = 1:k
    Jb = 0.5 * (J_blocks(:, :, ib) + J_blocks(:, :, ib).');
    Hb = 0.5 * (H_blocks(:, :, ib) + H_blocks(:, :, ib).');

    [Vb, lambdab, okb, Zb] = explicit_2x2_generalized_basis(Jb, Hb);
    methodb = 'partial-block-explicit';
    if ~okb && opt.safe_characteristic_fallback
        [Vb, lambdab, okb, Zb] = generalized_flux_basis(Jb, Hb);
        methodb = 'partial-block-generalized-fallback';
    end

    if ~okb
        ok = false;
        state = struct('method', methodb, 'block_methods', {methods});
        return;
    end

    idx = (2 * ib - 1):(2 * ib);
    V(idx, idx) = Vb;
    Z(idx, idx) = Zb;
    lambda(idx) = lambdab(:);
    methods{ib} = methodb;
end

state = struct('method', 'partial-block-explicit', 'block_methods', {methods});
if any(strcmp(methods, 'partial-block-generalized-fallback'))
    state.method = 'partial-block-mixed';
end
end

function [V, lambda, ok, Z] = explicit_2x2_generalized_basis(J, H)
V = eye(2);
Z = eye(2);
lambda = zeros(2, 1);
ok = false;

a = H(1, 1); b = H(1, 2); c = H(2, 2);
p = J(1, 1); q = J(1, 2); r = J(2, 2);

detH = a * c - b * b;
scale = max(1.0, norm(H, Inf) + norm(J, Inf));
tol = 1.0e-12 * scale;
if ~isfinite(detH) || detH <= tol
    return;
end

coeff2 = detH;
coeff1 = -p * c - a * r + 2.0 * q * b;
coeff0 = p * r - q * q;

disc = coeff1 * coeff1 - 4.0 * coeff2 * coeff0;
if disc < -tol * scale
    return;
end
disc = max(disc, 0.0);
sqrt_disc = sqrt(disc);

lambda = [(-coeff1 - sqrt_disc) / (2.0 * coeff2); ...
          (-coeff1 + sqrt_disc) / (2.0 * coeff2)];
if any(~isfinite(lambda))
    return;
end

for j = 1:2
    K = J - lambda(j) * H;
    row1 = K(1, :);
    row2 = K(2, :);
    if norm(row1, 2) >= norm(row2, 2)
        z = [-row1(2); row1(1)];
    else
        z = [-row2(2); row2(1)];
    end

    if norm(z, 2) <= tol
        z = [K(2, 2); -K(2, 1)];
    end
    if norm(z, 2) <= tol
        return;
    end

    Z(:, j) = z;
    V(:, j) = H * z;
end

if abs(det(V)) <= tol
    return;
end

ok = true;
end

function [V, lambda, ok, Z] = generalized_flux_basis(J, H)
n = size(H, 1);
V = eye(n);
Z = eye(n);
lambda = zeros(n, 1);
ok = false;

if any(~isfinite(J(:))) || any(~isfinite(H(:)))
    return;
end

[R, p] = chol(full(H));
if p ~= 0
    shift = max(1.0e-14, 1.0e-12 * max(1.0, norm(H, Inf)));
    [R, p] = chol(full(H + shift * speye(n)));
    if p ~= 0
        return;
    end
end

S = R' \ (full(J) / R);
S = 0.5 * (S + S.');

[Y, D] = eig(S);
lambda = diag(D);
if any(~isfinite(Y(:))) || any(~isfinite(lambda)) || any(abs(imag(lambda)) > 0)
    return;
end

Y = real(Y);
Z = R \ Y;
V = R' * Y;
ok = true;
end

function [V, lambda, ok] = nonsymmetric_real_basis(F)
n = size(F, 1);
V = zeros(n, n);
lambda = zeros(n, 1);
ok = false;

if any(~isfinite(F(:)))
    return;
end

[Z, D] = eig(F);
lam = diag(D);
if any(~isfinite(Z(:))) || any(~isfinite(lam))
    return;
end

lamScale = max(1.0, max(abs(real(lam))));
if any(abs(imag(lam)) > 1.0e-8 * lamScale)
    return;
end

[lamReal, perm] = sort(real(lam), 'ascend');
Z = Z(:, perm);

tol = 1.0e-8 * lamScale;
col = 1;
i = 1;
while i <= n
    j = i;
    while j < n && abs(lamReal(j + 1) - lamReal(i)) <= tol
        j = j + 1;
    end

    groupSize = j - i + 1;
    candidates = [real(Z(:, i:j)), imag(Z(:, i:j))];
    chosen = zeros(n, groupSize);
    nChosen = 0;
    for c = 1:size(candidates, 2)
        v = candidates(:, c);
        if norm(v) <= 1.0e-13
            continue;
        end
        for q = 1:nChosen
            v = v - chosen(:, q) * (chosen(:, q).' * v);
        end
        nv = norm(v);
        if nv > 1.0e-10
            nChosen = nChosen + 1;
            chosen(:, nChosen) = v / nv;
            if nChosen == groupSize
                break;
            end
        end
    end

    if nChosen < groupSize
        return;
    end

    V(:, col:(col + groupSize - 1)) = chosen;
    lambda(col:(col + groupSize - 1)) = lamReal(i:j);
    col = col + groupSize;
    i = j + 1;
end

ok = true;
end

function [V, Z, ok] = normalize_basis_columns(V, Z)
ok = false;
n = size(V, 2);
if nargin < 2
    Z = [];
end

for j = 1:n
    v = real(V(:, j));
    nv = norm(v);
    if ~isfinite(nv) || nv <= 0
        return;
    end
    v = v / nv;
    if ~isempty(Z)
        Z(:, j) = real(Z(:, j)) / nv;
    end

    [~, imax] = max(abs(v));
    if ~isempty(imax) && v(imax) < 0
        v = -v;
        if ~isempty(Z)
            Z(:, j) = -Z(:, j);
        end
    end
    V(:, j) = v;
end
ok = all(isfinite(V(:))) && (isempty(Z) || all(isfinite(Z(:))));
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

function residual = generalized_residual_norm(J, H, Z, lambda)
residual = norm(full(J * Z - H * (Z * diag(lambda))), Inf);
end

function residual = generalized_residual_block_norm(J_blocks, H_blocks, Z, lambda)
k = size(J_blocks, 3);
residual = 0.0;
for ib = 1:k
    idx = (2 * ib - 1):(2 * ib);
    Zb = Z(idx, idx);
    lambdab = lambda(idx);
    rb = J_blocks(:, :, ib) * Zb - H_blocks(:, :, ib) * (Zb * diag(lambdab));
    residual = max(residual, norm(rb, Inf));
end
end

function opt = characteristic_options(opt_cfg)
opt = struct();
opt.check_eigen_residual = logical(get_field_or(opt_cfg, 'check_eigen_residual', false));
opt.check_basis_condition = logical(get_field_or(opt_cfg, 'check_basis_condition', false));
opt.safe_characteristic_fallback = logical(get_field_or(opt_cfg, 'safe_characteristic_fallback', true));
opt.eigen_residual_tol = get_field_or(opt_cfg, 'eigen_residual_tol', 1.0e-7);
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
