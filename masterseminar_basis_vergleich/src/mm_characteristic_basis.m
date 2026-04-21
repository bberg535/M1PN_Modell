function [V, Vinv, state] = mm_characteristic_basis(u, model, quad_flux, opt_cfg, cache_state)
%MM_CHARACTERISTIC_BASIS Characteristic basis for flux reconstruction.
% Reference: Seminarquelle 2, Section 5.2.
% Computes the right eigenvectors of f'(u)=J(u)H(u)^{-1} through the
% symmetric generalized eigenproblem J ztilde = lambda H ztilde.

if nargin < 5
    cache_state = struct();
end

n = model.nMom;
V = eye(n);
Vinv = eye(n);
state = struct('success', false, 'method', 'identity', 'lambda', zeros(n, 1), ...
    'residual_inf', Inf);

try
    [J, H, matrix_state] = flux_jacobian_factors(u, model, quad_flux, opt_cfg, cache_state);
    state.matrix = matrix_state;

    if any(~isfinite(J(:))) || any(~isfinite(H(:)))
        return;
    end

    if strcmp(model.family, 'partial')
        [Vtry, lambda, ok] = partial_block_basis(J, H, model);
        method = 'partial-block-generalized';
    else
        [Vtry, lambda, ok] = generalized_flux_basis(J, H);
        method = 'symmetric-generalized';
    end

    if ~ok
        F = J / H;
        [Vtry, lambda, ok] = nonsymmetric_real_basis(F);
        method = 'nonsymmetric-real-fallback';
    end

    if ~ok
        return;
    end

    [lambda, idx] = sort(real(lambda(:)), 'ascend');
    Vtry = real(Vtry(:, idx));
    [Vtry, ok] = normalize_basis_columns(Vtry);
    if ~ok || rcond(Vtry) < 1.0e-12
        return;
    end

    VinvTry = Vtry \ eye(n);
    if any(~isfinite(VinvTry(:)))
        return;
    end

    residual = norm((J / H) * Vtry - Vtry * diag(lambda), Inf);
    scale = max(1.0, norm(J / H, Inf));
    if ~isfinite(residual) || residual > 1.0e-7 * scale
        return;
    end

    V = Vtry;
    Vinv = VinvTry;
    state.success = true;
    state.method = method;
    state.lambda = lambda;
    state.residual_inf = residual;
catch err
    state.error = err.message;
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
else
    weights = w;
    state.alpha = model.mass_matrix \ u(:);
    state.ansatz = struct('info', struct('converged', true, 'iterations', 0));
    state.closure = 'linear';
end

J = B.' * bsxfun(@times, weights .* mu, B);
if model.needs_entropy
    H = B.' * bsxfun(@times, weights, B);
else
    H = model.mass_matrix;
end

J = 0.5 * (J + J.');
H = 0.5 * (H + H.');
end

function [V, lambda, ok] = partial_block_basis(J, H, model)
n = model.nMom;
V = zeros(n, n);
lambda = zeros(n, 1);
ok = true;

for ib = 1:model.kIntervals
    idx = (2 * ib - 1):(2 * ib);
    [Vb, lambdab, okb] = generalized_flux_basis(J(idx, idx), H(idx, idx));
    if ~okb
        ok = false;
        return;
    end
    V(idx, idx) = Vb;
    lambda(idx) = lambdab;
end
end

function [V, lambda, ok] = generalized_flux_basis(J, H)
n = size(J, 1);
V = eye(n);
lambda = zeros(n, 1);
ok = false;

if any(~isfinite(J(:))) || any(~isfinite(H(:)))
    return;
end

[R, p] = chol(H);
if p ~= 0
    shift = max(1.0e-14, 1.0e-12 * max(1.0, norm(H, Inf)));
    [R, p] = chol(H + shift * eye(n));
    if p ~= 0
        return;
    end
end

S = R' \ (J / R);
S = 0.5 * (S + S.');

[Y, D] = eig(S);
lambda = diag(D);
if any(~isfinite(Y(:))) || any(~isfinite(lambda)) || any(abs(imag(lambda)) > 0)
    return;
end

% If ztilde solves J ztilde=lambda H ztilde, then z=H ztilde=R'Y is a
% right eigenvector of JH^{-1}.
V = R' * real(Y);
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

function [V, ok] = normalize_basis_columns(V)
ok = false;
n = size(V, 2);
for j = 1:n
    v = real(V(:, j));
    nv = norm(v);
    if ~isfinite(nv) || nv <= 0
        return;
    end
    v = v / nv;
    [~, imax] = max(abs(v));
    if ~isempty(imax) && v(imax) < 0
        v = -v;
    end
    V(:, j) = v;
end
ok = all(isfinite(V(:)));
end
