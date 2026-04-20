function ref = mm_reference_sn_source_beam(cfg_ref)
%MM_REFERENCE_SN_SOURCE_BEAM High-resolution S_N reference for source-beam benchmark.
% Benchmark setup follows Seminarquelle 2, Section 6.1.2.

if nargin < 1
    cfg_ref = struct();
end

domain = get_field_or(cfg_ref, 'domain', [0.0, 3.0]);
zL = domain(1);
zR = domain(2);
tf = get_field_or(cfg_ref, 'tf', 2.5);
nCells = get_field_or(cfg_ref, 'n_cells', 1200);
nMu = get_field_or(cfg_ref, 'n_mu', 128);
cfl = get_field_or(cfg_ref, 'cfl', 0.45);
psi_vac = get_field_or(cfg_ref, 'psi_vac_density', 1.0e-8 / 2.0);
beamExp = get_field_or(cfg_ref, 'beam_exponent', 1.0e5);
normalizeLeft = logical(get_field_or(cfg_ref, 'left_boundary_normalize', false));
par_cfg = get_field_or(cfg_ref, 'parallel', struct());

edges = linspace(zL, zR, nCells + 1);
z = 0.5 * (edges(1:end-1) + edges(2:end));
dz = edges(2) - edges(1);

sigma_a = eval_cell_field(get_field_or(cfg_ref, 'sigma_a', []), z, @default_sigma_a);
sigma_s = eval_cell_field(get_field_or(cfg_ref, 'sigma_s', []), z, @default_sigma_s);
Q = eval_cell_field(get_field_or(cfg_ref, 'Q', []), z, @default_Q);

[mu, w] = gauss_legendre(nMu, -1, 1);
maxSpeed = max(abs(mu));
par_mu = mm_parallel_context(par_cfg, nMu, 'directions');
use_par_mu = par_mu.use_parallel;

left_fun = get_field_or(cfg_ref, 'boundary_psi_left', ...
    @(mu_vals, t) default_left_bc(mu_vals, t, beamExp, normalizeLeft)); %#ok<NASGU>
right_fun = get_field_or(cfg_ref, 'boundary_psi_right', ...
    @(mu_vals, t) psi_vac * ones(size(mu_vals))); %#ok<NASGU>

psi = psi_vac * ones(nMu, nCells);

t = 0.0;
while t < tf - 1e-14
    dt = min(cfl * dz / maxSpeed, tf - t);

    k1 = rhs_sn(psi, t);
    psi1 = max(psi + dt * k1, 0.0);

    k2 = rhs_sn(psi1, t + dt);
    psi = max(psi + 0.5 * dt * (k1 + k2), 0.0);

    t = t + dt;
end

rho = w.' * psi;

ref = struct();
ref.z = z(:);
ref.rho = rho(:);
ref.t = t;
ref.mu = mu;
ref.w = w;

    function rhs = rhs_sn(psi_state, t_state)
        rho_state = w.' * psi_state;
        src = bsxfun(@times, sigma_s, 0.5 * rho_state - psi_state) - ...
              bsxfun(@times, sigma_a, psi_state) + ...
              repmat(Q, nMu, 1);

        psi_in_left = eval_boundary_fun(left_fun, mu, t_state, psi_vac);
        psi_in_right = eval_boundary_fun(right_fun, mu, t_state, psi_vac);

        adv = zeros(size(psi_state));
        if use_par_mu
            parfor m = 1:nMu
                mmu = mu(m);
                row = psi_state(m, :);

                if mmu >= 0
                    F = zeros(1, nCells + 1);
                    F(1) = mmu * psi_in_left(m);
                    F(2:end) = mmu * row;
                else
                    F = zeros(1, nCells + 1);
                    F(1:nCells) = mmu * row;
                    F(end) = mmu * psi_in_right(m);
                end

                adv(m, :) = -(F(2:end) - F(1:end-1)) / dz;
            end
        else
            for m = 1:nMu
                mmu = mu(m);
                row = psi_state(m, :);

                if mmu >= 0
                    F = zeros(1, nCells + 1);
                    F(1) = mmu * psi_in_left(m);
                    F(2:end) = mmu * row;
                else
                    F = zeros(1, nCells + 1);
                    F(1:nCells) = mmu * row;
                    F(end) = mmu * psi_in_right(m);
                end

                adv(m, :) = -(F(2:end) - F(1:end-1)) / dz;
            end
        end

        rhs = adv + src;
    end

end

function v = eval_cell_field(spec, z, default_fun)
if isempty(spec)
    v = default_fun(z);
elseif isscalar(spec)
    v = spec * ones(size(z));
elseif isa(spec, 'function_handle')
    v = spec(z);
else
    v = spec;
end
v = reshape(v, 1, []);
if numel(v) ~= numel(z)
    error('Source-beam field has invalid size.');
end
end

function sigma_a = default_sigma_a(z)
sigma_a = zeros(size(z));
sigma_a(z <= 2.0) = 1.0;
end

function sigma_s = default_sigma_s(z)
sigma_s = zeros(size(z));
sigma_s(z > 1.0 & z <= 2.0) = 2.0;
sigma_s(z > 2.0) = 10.0;
end

function q = default_Q(z)
q = zeros(size(z));
q(z >= 1.0 & z <= 1.5) = 0.5;
end

function psi = default_left_bc(mu, ~, beamExp, normalize_flag)
psi = exp(-beamExp * (mu - 1.0).^2);
if normalize_flag
    denom = trapz(mu(mu >= 0), psi(mu >= 0));
    if denom > 0
        psi = psi / denom;
    end
end
end

function y = eval_boundary_fun(fun, mu, t, fallback)
if isa(fun, 'function_handle')
    try
        y = fun(mu, t);
    catch
        try
            y = fun(mu);
        catch
            try
                y = fun(t, mu);
            catch
                y = fun(t);
            end
        end
    end
else
    y = fun;
end

if isscalar(y)
    y = y * ones(size(mu));
end
y = reshape(y, size(mu));
y = max(y, 0.0);

if nargin >= 4 && isempty(y)
    y = fallback * ones(size(mu));
end
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
end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
