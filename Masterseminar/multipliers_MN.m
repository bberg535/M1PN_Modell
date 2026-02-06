function [alpha_cache, ok] = multipliers_MN(u, alpha_cache, opts, Pfull)
% Solve u = < b exp(alpha·b) > cellwise using Newton.
% Uses density scaling idea (paper) but implemented in a compact, robust way.

nM = size(u,1);
Nx = size(u,2);
quad = opts.quad;
w = quad.w_full(:);
ok = true;

for i=1:Nx
    ui = u(:,i);
    rho = ui(1);

    % enforce minimal density (regularization idea in the paper)
    if rho < opts.rho_vac
        ui = zeros(nM,1); ui(1)=opts.rho_vac;
        u(:,i) = ui;
        alpha_cache(:,i) = [log(opts.rho_vac/2); zeros(nM-1,1)];
        continue;
    end

    % scale to density ~1 for Newton (paper-style)
    phi = ui / rho;

    % initial guess for beta (density~1)
    beta = alpha_cache(:,i);
    % shift guess to roughly match scaled density:
    beta(1) = beta(1) - log(max(rho, 1e-300));

    % Newton on F(beta)=u(beta)-phi
    [beta, success] = newton_beta(phi, beta, Pfull, w, opts.newton_tol, opts.newton_maxit, opts.newton_damp);
    if ~success
        ok = false;
        % fallback: isotropic
        beta = zeros(nM,1);
        beta(1) = log(1/2);
    end

    % rescale back to preserve exact rho (paper formula concept)
    ubeta = moment_map(beta, Pfull, w);
    rho_beta = max(ubeta(1), 1e-300);

    alpha = beta;
    alpha(1) = alpha(1) + log(rho/rho_beta);

    alpha_cache(:,i) = alpha;
end
end
