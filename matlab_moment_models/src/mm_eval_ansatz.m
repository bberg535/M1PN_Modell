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
    [alpha, info] = mm_entropy_dual_solve(u, model, quad_eval, opt_cfg, cache_state);
    psi = exp(min(B * alpha, 700));
    aux = struct('info', info);
else
    alpha = model.mass_matrix \ u;
    psi = B * alpha;
    aux = struct('info', struct('converged', true, 'iterations', 0));
end

end
