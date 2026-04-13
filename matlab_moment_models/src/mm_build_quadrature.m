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
if model.is_partial
    quad.B = build_partial_basis_with_interval_ids(model, quad.mu, interval_id(:));
    quad.interval_index = interval_id(:);
else
    quad.B = model.basis_eval(quad.mu);
end

pos = quad.mu >= 0;
neg = quad.mu <= 0;
quad.mu_plus = quad.mu(pos);
quad.w_plus = quad.w(pos);
quad.B_plus = quad.B(pos, :);
quad.mu_minus = quad.mu(neg);
quad.w_minus = quad.w(neg);
quad.B_minus = quad.B(neg, :);

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

function B = build_partial_basis_with_interval_ids(model, mu, interval_id)
nQ = numel(mu);
B = zeros(nQ, model.nMom);

for q = 1:nQ
    j = interval_id(q);
    if j < 1 || j > model.kIntervals
        continue;
    end
    col = 2 * j - 1;
    B(q, col) = 1.0;
    B(q, col + 1) = mu(q);
end
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
