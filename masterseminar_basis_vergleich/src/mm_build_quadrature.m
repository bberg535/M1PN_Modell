function quad = mm_build_quadrature(model, cfg_quad, purpose)
%MM_BUILD_QUADRATURE Build quadrature for moments/flux/realizability checks.
% References:
% - Quadrature strategy in slab: Seminarquelle 2, Section 5.4
% - Half-interval split [-1,0]/[0,1] for kinetic flux integrals in Eq. (4.15)

if nargin < 3 || isempty(purpose)
    purpose = 'default';
end

purpose = lower(strtrim(purpose));
quad = struct();
quad.purpose = purpose;

switch model.family
    case {'hat', 'partial'}
        % Gauss-Lobatto per interval (Section 5.4, suitable for HFM/PMM realizability).
        qOrder = get_field_or(cfg_quad, 'lobatto_order', 15);
        [mu, w, interval_id] = piecewise_lobatto(model.mu_edges, qOrder);

    case 'full'
        % For full moments (PN/MN), use high-order quadrature on both half intervals.
        baseOrder = get_field_or(cfg_quad, 'base_legendre_order', 80);
        fluxOrder = max(baseOrder, 2 * model.N + 40);
        switch purpose
            case {'flux', 'lp', 'default', 'moment'}
                [muL, wL] = gauss_legendre(fluxOrder, -1, 0);
                [muR, wR] = gauss_legendre(fluxOrder, 0, 1);
                mu = [muL; muR];
                w = [wL; wR];
            otherwise
                [mu, w] = gauss_legendre(fluxOrder, -1, 1);
        end

    otherwise
        error('Unknown model family: %s', model.family);
end

quad.mu = mu(:);
quad.w = w(:);
if ismember(model.family, {'hat', 'partial'})
    quad.interval_id = interval_id(:);
    quad.interval_index = quad.interval_id;
    [quad.B, quad.local_mu_by_interval, quad.local_w_by_interval, ...
        quad.local_B_by_interval, quad.local_cols_by_interval, quad.interval_qidx] = ...
        build_local_piecewise_basis(model, quad.mu, quad.w, quad.interval_id);
else
    quad.B = model.basis_eval(quad.mu);
end
quad.Bt_full = full(quad.B.');
if strcmp(model.family, 'full') && strcmp(model.realizability, 'lp')
    quad.lp_hull = build_fullmoment_lp_hull(quad.Bt_full);
end

pos = quad.mu >= 0;
neg = quad.mu <= 0;
quad.mu_plus = quad.mu(pos);
quad.w_plus = quad.w(pos);
quad.B_plus = quad.B(pos, :);
quad.mu_minus = quad.mu(neg);
quad.w_minus = quad.w(neg);
quad.B_minus = quad.B(neg, :);

if ~model.needs_entropy
    quad.linear_char = build_linear_characteristic_data(model, quad);
end

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end

function [mu, w, interval_id] = piecewise_lobatto(edges, qOrder)
mu = [];
w = [];
interval_id = [];
for j = 1:(numel(edges) - 1)
    [mj, wj] = gauss_lobatto(qOrder, edges(j), edges(j + 1));
    mu = [mu; mj(:)]; %#ok<AGROW>
    w = [w; wj(:)]; %#ok<AGROW>
    interval_id = [interval_id; j * ones(numel(mj), 1)]; %#ok<AGROW>
end
end

function [B, local_mu, local_w, local_B, local_cols, interval_qidx] = build_local_piecewise_basis(model, mu, w, interval_id)
nQ = numel(mu);
nIntervals = numel(model.mu_edges) - 1;
row_idx = [];
col_idx = [];
val = [];

local_mu = cell(1, nIntervals);
local_w = cell(1, nIntervals);
local_B = cell(1, nIntervals);
local_cols = cell(1, nIntervals);
interval_qidx = cell(1, nIntervals);

for j = 1:nIntervals
    qidx = find(interval_id == j);
    interval_qidx{j} = qidx(:);
    mu_j = mu(qidx);
    w_j = w(qidx);
    local_mu{j} = mu_j(:);
    local_w{j} = w_j(:);

    switch model.family
        case 'partial'
            cols = [2 * j - 1, 2 * j];
            Bj = [ones(numel(mu_j), 1), mu_j(:)];

        case 'hat'
            cols = [j, j + 1];
            a = model.mu_edges(j);
            b = model.mu_edges(j + 1);
            h = b - a;
            Bj = [(b - mu_j(:)) / h, (mu_j(:) - a) / h];

        otherwise
            error('Unsupported piecewise family: %s', model.family);
    end

    local_B{j} = Bj;
    local_cols{j} = cols;

    if ~isempty(qidx)
        nLoc = numel(qidx);
        col_block = repmat(cols(:).', nLoc, 1);
        row_block = repmat(qidx(:), 1, numel(cols));
        row_idx = [row_idx; row_block(:)]; %#ok<AGROW>
        col_idx = [col_idx; col_block(:)]; %#ok<AGROW>
        val = [val; Bj(:)]; %#ok<AGROW>
    end
end

B = sparse(row_idx, col_idx, val, nQ, model.nMom);
end

function char_data = build_linear_characteristic_data(model, quad)
char_data = struct('constant', true, 'success', false, 'method', 'linear-generalized', ...
    'V', [], 'Vinv', [], 'lambda', [], 'Z', [], 'K_flux', [], 'M', model.mass_matrix);

K_flux = weighted_basis_gram(quad.B, quad.w(:) .* quad.mu(:));
[V, Vinv, lambda, Z, ok] = solve_generalized_constant(K_flux, model.mass_matrix);

char_data.K_flux = K_flux;
char_data.J = K_flux;
char_data.H = model.mass_matrix;
if ok
    char_data.success = true;
    char_data.V = V;
    char_data.Vinv = Vinv;
    char_data.lambda = lambda;
    char_data.Z = Z;
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

function [V, Vinv, lambda, Z, ok] = solve_generalized_constant(J, H)
n = size(H, 1);
V = eye(n);
Vinv = eye(n);
lambda = zeros(n, 1);
Z = eye(n);
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
lambda = real(diag(D));
if any(~isfinite(Y(:))) || any(~isfinite(lambda))
    return;
end

[lambda, idx] = sort(lambda, 'ascend');
Z = R \ real(Y(:, idx));
V = real(R' * Y(:, idx));
[V, Z, ok] = normalize_basis_columns(V, Z);
if ~ok
    return;
end

Vinv = V \ eye(n);
ok = all(isfinite(Vinv(:)));
end

function [V, Z, ok] = normalize_basis_columns(V, Z)
ok = false;
for j = 1:size(V, 2)
    v = real(V(:, j));
    nv = norm(v);
    if ~isfinite(nv) || nv <= 0
        return;
    end
    v = v / nv;
    Z(:, j) = Z(:, j) / nv;
    [~, imax] = max(abs(v));
    if ~isempty(imax) && v(imax) < 0
        v = -v;
        Z(:, j) = -Z(:, j);
    end
    V(:, j) = v;
end
ok = all(isfinite(V(:))) && all(isfinite(Z(:)));
end

function hull = build_fullmoment_lp_hull(Bt_full)
points = Bt_full(2:end, :).';
points = unique(points, 'rows', 'stable');
d = size(points, 2);
hull = struct('dimension', d, 'normals', [], 'offsets', [], 'bounds', []);

if d == 0
    return;
end

if d == 1
    hull.bounds = [min(points(:, 1)), max(points(:, 1))];
    return;
end

facets = convhulln(points);
center = mean(points, 1);
nFacets = size(facets, 1);
normals = zeros(nFacets, d);
offsets = zeros(nFacets, 1);

for i = 1:nFacets
    verts = points(facets(i, :), :);
    base = verts(1, :).';
    D = (verts(2:end, :) - verts(1, :)).';
    [~, ~, V] = svd(D.', 'econ');
    nvec = V(:, end);
    nvec = nvec / norm(nvec);
    offset = nvec.' * base;

    if nvec.' * center(:) > offset
        nvec = -nvec;
        offset = -offset;
    end

    normals(i, :) = nvec(:).';
    offsets(i) = offset;
end

hull.normals = normals;
hull.offsets = offsets;
end

function [x, w] = gauss_lobatto(n, a, b)
if n < 2
    error('Gauss-Lobatto requires at least n=2 points.');
end
if n == 2
    x = [a; b];
    w = (b - a) * [0.5; 0.5];
    return;
end

N = n - 1;
k = (0:N)';
x = cos(pi * k / N);

xold = 2 * ones(size(x));
P = zeros(n, n);
iter = 0;
while max(abs(x - xold)) > 1e-14 && iter < 200
    xold = x;
    P(:, 1) = 1.0;
    P(:, 2) = x;
    for m = 2:N
        P(:, m + 1) = ((2*m - 1) .* x .* P(:, m) - (m - 1) .* P(:, m - 1)) / m;
    end
    x = xold - (x .* P(:, N + 1) - P(:, N)) ./ (N * (N + 1) * P(:, N + 1));
    iter = iter + 1;
end

w = 2.0 ./ (N * (N + 1) * (P(:, N + 1).^2));

% map from [-1,1] to [a,b]
x = 0.5 * ((b - a) .* x + (a + b));
w = 0.5 * (b - a) .* w;

[x, idx] = sort(x);
w = w(idx);
end

function [x, w] = gauss_legendre(n, a, b)
i = (1:n-1)';
beta = i ./ sqrt(4 * i.^2 - 1);
T = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(T);
x = diag(D);
[x, idx] = sort(x);
V = V(:, idx);
w = 2 * (V(1, :)').^2;

x = 0.5 * ((b - a) * x + (a + b));
w = 0.5 * (b - a) * w;

x = x(:);
w = w(:);
end
