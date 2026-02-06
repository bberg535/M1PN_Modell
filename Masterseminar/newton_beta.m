function [beta, success] = newton_beta(phi, beta, Pfull, w, tol, maxit, damp)
success = true;

for it=1:maxit
    ub = moment_map(beta, Pfull, w);
    F  = ub - phi;

    if norm(F,2) < tol
        return;
    end

    H = hessian_map(beta, Pfull, w);

    % Newton direction
    d = H \ F;

    % damping / backtracking on residual norm
    t = 1.0;
    nF = norm(F,2);
    for ls=1:20
        beta_try = beta - t*d;
        ub_try = moment_map(beta_try, Pfull, w);
        F_try = ub_try - phi;
        if isfinite(norm(F_try,2)) && norm(F_try,2) < (1 - 1e-4*t)*nF
            beta = beta_try;
            break;
        end
        t = t*damp;
    end

    if t < 1e-12
        success = false;
        return;
    end
end

success = false;
end

