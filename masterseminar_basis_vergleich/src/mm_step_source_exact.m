function [u_next, source_state] = mm_step_source_exact(u, model, phys, dt, source_cfg)
%MM_STEP_SOURCE_EXACT Exact source substep for isotropic-scattering slab case.
% References (Seminarquelle 2):
% - Exact source evolution formula: Eq. (4.11)
% - Isotropic-source simplification used here: Eq. (4.12)

if nargin < 5
    source_cfg = struct();
end
tTotal = tic;

[nMom, nCells] = size(u);

sigma_s = expand_to_cells(phys.sigma_s, nCells);
sigma_a = expand_to_cells(phys.sigma_a, nCells);

if isfield(phys, 'Q')
    Q = expand_to_cells(phys.Q, nCells);
else
    Q = zeros(1, nCells);
end

rho = model.alpha1.' * u;
u_iso = model.u_iso_vec * rho;

exp_ss = exp(-sigma_s * dt);
exp_sa = exp(-sigma_a * dt);

% Homogeneous part of Eq. (4.12).
term_hom = exp_sa .* (exp_ss .* u + (1 - exp_ss) .* u_iso);

hbQ = model.b_iso * Q;
source_factor = zeros(1, nCells);
mask = abs(sigma_a) > 1e-14;
source_factor(mask) = (1 - exp_sa(mask)) ./ sigma_a(mask);
source_factor(~mask) = dt;

% Source contribution of Eq. (4.12), including sigma_a -> 0 limit.
term_src = hbQ .* source_factor;

u_next = term_hom + term_src;

source_state = struct();
source_state.rho = rho;
source_state.u_iso = u_iso;
source_state.source_factor = source_factor;

if isfield(source_cfg, 'enforce_realizability') && source_cfg.enforce_realizability
    for i = 1:nCells
        [ok, ~] = mm_is_realizable(u_next(:, i), model, mm_build_quadrature(model, struct(), 'lp'), struct('epsR', 0));
        if ~ok
            u_next(:, i) = u_iso(:, i);
        end
    end
end

source_state.timing = struct();
source_state.timing.total_s = toc(tTotal);

end

function arr = expand_to_cells(v, nCells)
if isscalar(v)
    arr = repmat(v, 1, nCells);
else
    arr = reshape(v, 1, []);
    if numel(arr) ~= nCells
        error('Field size does not match number of cells.');
    end
end
end
