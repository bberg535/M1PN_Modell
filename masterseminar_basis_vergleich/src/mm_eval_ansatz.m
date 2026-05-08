function [psi, alpha, aux] = mm_eval_ansatz(u, model, quad_eval, opt_cfg, cache_state)
%MM_EVAL_ANSATZ Evaluate closure ansatz density at quadrature nodes.
% References:
% - Entropy ansatz reconstruction: Seminarquelle 1 Eq. (2.6), Seminarquelle 2 Eq. (2.5)/(2.6).
% - Linear ansatz psi=b^T a(u): Seminarquelle 1 Eq. (2.12).

if nargin < 5
    cache_state = struct();
end

u = u(:);
B = quad_eval.B;

if model.needs_entropy
    if isfield(cache_state, 'psi') && numel(cache_state.psi) == size(B, 1)
        psi = cache_state.psi(:);
        alpha = get_field_or(cache_state, 'alpha', zeros(model.nMom, 1));
        aux = struct('info', get_field_or(cache_state, 'info', struct('converged', true, 'iterations', 0)), ...
            'from_cache', true);
        return;
    end

    if isfield(cache_state, 'alpha') && numel(cache_state.alpha) == model.nMom
        alpha = cache_state.alpha(:);
        psi = exp(min(B * alpha, 700));
        aux = struct('info', get_field_or(cache_state, 'info', struct('converged', true, 'iterations', 0)), ...
            'from_cache', true);
        return;
    end

    [alpha, info] = mm_entropy_dual_solve(u, model, quad_eval, opt_cfg, cache_state);
    psi = exp(min(B * alpha, 700));
    aux = struct('info', info, 'from_cache', false);
else
    if isfield(cache_state, 'alpha') && numel(cache_state.alpha) == model.nMom
        alpha = cache_state.alpha(:);
    else
        alpha = model.mass_matrix \ u;
    end
    if isfield(cache_state, 'psi') && numel(cache_state.psi) == size(B, 1)
        psi = cache_state.psi(:);
        from_cache = true;
    else
        psi = B * alpha;
        from_cache = isfield(cache_state, 'alpha');
    end
    aux = struct('info', struct('converged', true, 'iterations', 0), 'from_cache', from_cache);
end

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
