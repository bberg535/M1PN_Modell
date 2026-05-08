function [f, state] = mm_eval_flux_function(u, model, quad_flux, opt_cfg, cache_state)
%MM_EVAL_FLUX_FUNCTION Evaluate physical flux f(u)=<mu b psi_u>.
% Reference: Seminarquelle 2 Eq. (2.13) (slab moment system flux term).

if nargin < 5
    cache_state = struct();
end

if isfield(cache_state, 'flux') && numel(cache_state.flux) == model.nMom
    f = cache_state.flux(:);
    state = struct('alpha', get_field_or(cache_state, 'alpha', []), ...
        'aux', struct('from_cache', true, 'info', get_field_or(cache_state, 'info', struct())), ...
        'psi', get_field_or(cache_state, 'psi', []), ...
        'block_moments', get_field_or(cache_state, 'block_moments', []), ...
        'from_cache', true);
    return;
end

[psi, alpha, aux] = mm_eval_ansatz(u, model, quad_flux, opt_cfg, cache_state);

f = quad_flux.B.' * (quad_flux.w .* quad_flux.mu .* psi);
state = struct('alpha', alpha, 'aux', aux, 'psi', psi, 'block_moments', [], 'from_cache', false);

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
