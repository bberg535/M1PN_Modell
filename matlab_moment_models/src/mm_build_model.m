function model = mm_build_model(model_name, order, cfg_model)
%MM_BUILD_MODEL Build 1D slab angular model descriptor.
% References (Seminarquelle 1):
% - Full-moment bases f_N: Eq. (3.1)/(3.2)
% - Hat basis h_i: Eq. (3.3a)
% - Linear vs entropy closure: Eq. (2.11)/(2.10), ansatz Eq. (2.12)

if nargin < 3
    cfg_model = struct();
end

name = normalize_name(model_name);

model = struct();
model.name = name;
model.order = order;
model.ndim = 1;
model.domain = [-1, 1];
model.cfg = cfg_model;

switch name
    case {'PN', 'MN'}
        % We use the Legendre full-moment basis f_N (Eq. (3.2)).
        N = order;
        n = N + 1;
        model.family = 'full';
        model.nMom = n;
        model.N = N;
        model.mu_edges = [];
        model.alpha1 = zeros(n, 1);
        model.alpha1(1) = 1.0;
        model.b_iso = zeros(n, 1);
        model.b_iso(1) = 2.0;
        model.mass_matrix = diag(2.0 ./ (2.0 * (0:N)' + 1.0));
        model.realizability = 'lp';
        model.is_partial = false;
        model.is_hat = false;

    case {'HFPN', 'HFMN'}
        % Piecewise linear hat basis h_i on interval partition, Eq. (3.3a).
        n = order;
        if n < 2
            error('Hat basis requires n >= 2.');
        end
        mu_edges = linspace(-1, 1, n);
        model.family = 'hat';
        model.nMom = n;
        model.N = n - 1;
        model.mu_edges = mu_edges(:);
        model.alpha1 = ones(n, 1);
        model.b_iso = hat_integrals(mu_edges(:));
        model.mass_matrix = hat_mass_matrix(mu_edges(:));
        model.realizability = 'positive';
        model.is_partial = false;
        model.is_hat = true;

    case {'PMPN', 'PMMN'}
        % Piecewise first-order partial moments p_{I_j}, Definition 3.6.
        n = order;
        if mod(n, 2) ~= 0 || n < 2
            error('Partial moment basis requires even n >= 2.');
        end
        k = n / 2;
        mu_edges = linspace(-1, 1, k + 1);
        model.family = 'partial';
        model.nMom = n;
        model.N = 1;
        model.kIntervals = k;
        model.mu_edges = mu_edges(:);
        model.alpha1 = zeros(n, 1);
        model.alpha1(1:2:end) = 1.0;
        model.b_iso = partial_integrals(mu_edges(:));
        model.mass_matrix = partial_mass_matrix(mu_edges(:));
        model.realizability = 'interval';
        model.is_partial = true;
        model.is_hat = false;

    otherwise
        error('Unknown model name: %s', model_name);
end

model.h1 = 2.0;
model.u_iso_vec = model.b_iso / model.h1;
model.G = model.u_iso_vec * model.alpha1.';
model.needs_entropy = ismember(name, {'MN', 'HFMN', 'PMMN'});

if model.needs_entropy
    % Entropy closure with Maxwell-Boltzmann entropy, Eq. (2.10).
    model.closure = 'entropy';
else
    % Linear closure using quadratic entropy, Eq. (2.11) and ansatz Eq. (2.12).
    model.closure = 'linear';
end

if ismember(name, {'MN', 'PMMN'})
    model.use_change_of_basis_default = true;
elseif strcmp(name, 'HFMN')
    model.use_change_of_basis_default = false;
else
    model.use_change_of_basis_default = false;
end

model.basis_eval = @(mu) basis_eval(mu, model);

end

function name = normalize_name(s)
s = upper(strtrim(s));
s = strrep(s, ' ', '');
if strcmp(s, 'HFP')
    s = 'HFPN';
elseif strcmp(s, 'HFM')
    s = 'HFMN';
elseif strcmp(s, 'PMP')
    s = 'PMPN';
elseif strcmp(s, 'PMM')
    s = 'PMMN';
end
name = s;
end

function B = basis_eval(mu, model)
mu = mu(:);
B = zeros(numel(mu), model.nMom);

switch model.family
    case 'full'
        N = model.N;
        B(:, 1) = 1.0;
        if N >= 1
            B(:, 2) = mu;
        end
        for l = 2:N
            % Legendre recursion: P_l(x)
            B(:, l + 1) = ((2*l - 1) .* mu .* B(:, l) - (l - 1) .* B(:, l - 1)) / l;
        end

    case 'hat'
        % Hat basis evaluation according to Eq. (3.3a).
        e = model.mu_edges;
        n = model.nMom;
        for i = 1:n
            if i == 1
                a = e(1); b = e(2);
                idx = mu >= a & mu <= b;
                B(idx, i) = (b - mu(idx)) / (b - a);
            elseif i == n
                a = e(end - 1); b = e(end);
                idx = mu >= a & mu <= b;
                B(idx, i) = (mu(idx) - a) / (b - a);
            else
                a = e(i - 1); b = e(i); c = e(i + 1);
                idx1 = mu >= a & mu <= b;
                idx2 = mu >= b & mu <= c;
                B(idx1, i) = (mu(idx1) - a) / (b - a);
                B(idx2, i) = (c - mu(idx2)) / (c - b);
            end
        end

    case 'partial'
        % Local partial-moment basis p_{I_j}=(1,mu) on each interval I_j.
        e = model.mu_edges;
        k = model.kIntervals;
        for j = 1:k
            a = e(j); b = e(j + 1);
            idx = mu >= a & mu <= b;
            col = 2*j - 1;
            B(idx, col) = 1.0;
            B(idx, col + 1) = mu(idx);
        end

    otherwise
        error('Unsupported model family: %s', model.family);
end

end

function b_iso = hat_integrals(edges)
n = numel(edges);
b_iso = zeros(n, 1);
for i = 1:n
    if i == 1
        h = edges(2) - edges(1);
        b_iso(i) = h / 2.0;
    elseif i == n
        h = edges(end) - edges(end - 1);
        b_iso(i) = h / 2.0;
    else
        hL = edges(i) - edges(i - 1);
        hR = edges(i + 1) - edges(i);
        b_iso(i) = (hL + hR) / 2.0;
    end
end
end

function M = hat_mass_matrix(edges)
n = numel(edges);
M = zeros(n, n);
for i = 1:(n - 1)
    h = edges(i + 1) - edges(i);
    M(i, i) = M(i, i) + h / 3.0;
    M(i + 1, i + 1) = M(i + 1, i + 1) + h / 3.0;
    M(i, i + 1) = M(i, i + 1) + h / 6.0;
    M(i + 1, i) = M(i + 1, i) + h / 6.0;
end
end

function b_iso = partial_integrals(edges)
k = numel(edges) - 1;
b_iso = zeros(2*k, 1);
for j = 1:k
    a = edges(j);
    b = edges(j + 1);
    idx = 2*j - 1;
    b_iso(idx) = b - a;
    b_iso(idx + 1) = 0.5 * (b^2 - a^2);
end
end

function M = partial_mass_matrix(edges)
k = numel(edges) - 1;
M = zeros(2*k, 2*k);
for j = 1:k
    a = edges(j);
    b = edges(j + 1);
    idx = 2*j - 1;
    m00 = b - a;
    m01 = 0.5 * (b^2 - a^2);
    m11 = (b^3 - a^3) / 3.0;
    M(idx:idx+1, idx:idx+1) = [m00, m01; m01, m11];
end
end
