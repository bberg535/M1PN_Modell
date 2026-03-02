function [f, state] = mm_eval_flux_function(u, model, quad_flux, opt_cfg, cache_state)
%MM_EVAL_FLUX_FUNCTION Evaluate physical flux f(u)=<mu b psi_u>.
% Reference: Seminarquelle 2 Eq. (2.13) (slab moment system flux term).

if nargin < 5
    cache_state = struct();
end

[psi, alpha, aux] = mm_eval_ansatz(u, model, quad_flux, opt_cfg, cache_state);

f = quad_flux.B.' * (quad_flux.w .* quad_flux.mu .* psi);
state = struct('alpha', alpha, 'aux', aux);

end
