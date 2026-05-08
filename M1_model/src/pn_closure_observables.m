function obs = pn_closure_observables(uPN, uM1)
%PN_CLOSURE_OBSERVABLES Compute closure-consistent PN observables for current M1 moments.

    nec = size(uM1, 2);
    [xi, w, ~, L] = pn_basis(size(uPN, 1) - 1);
    A = [w.'; (w(:) .* xi(:)).'];

    Vtrial = L * uPN(:, 1:nec);
    B = uM1(:, 1:nec);
    Vfit = Vtrial - A.' * ((A * A.') \ (A * Vtrial - B));

    [Vref, feasible] = nonnegative_reference(B, xi, w);
    if any(feasible)
        theta = ones(1, nec);
        for k = 1:size(Vfit, 1)
            mask = feasible & (Vfit(k, :) < 0);
            if any(mask)
                denom = Vref(k, mask) - Vfit(k, mask);
                theta(mask) = min(theta(mask), Vref(k, mask) ./ max(denom, eps));
            end
        end
        theta = max(0, min(1, theta));
        Vfit(:, feasible) = Vref(:, feasible) ...
            + (Vfit(:, feasible) - Vref(:, feasible)) .* theta(feasible);
    end

    obs = struct();
    obs.rho = B(1, :).';
    obs.j = B(2, :).';
    obs.m2 = ((w(:) .* (xi(:) .^ 2)).' * Vfit).';
    obs.angular_values = Vfit;
    obs.min_reconstruction = min(Vfit, [], 1).';
end

function [Vref, feasible] = nonnegative_reference(B, xi, w)
    rho = B(1, :);
    mom = B(2, :);
    K = numel(xi);
    nec = size(B, 2);
    tol = 1e-12;

    Vref = zeros(K, nec);
    feasible = rho <= tol;
    active = rho > tol;
    reduced_flux = zeros(1, nec);
    reduced_flux(active) = mom(active) ./ rho(active);

    assigned = false(1, nec);
    for j = 1:K-1
        mask = active & ~assigned ...
            & (reduced_flux >= xi(j) - tol) ...
            & (reduced_flux <= xi(j + 1) + tol);
        if any(mask)
            denom = xi(j + 1) - xi(j);
            Vref(j, mask) = rho(mask) .* (xi(j + 1) - reduced_flux(mask)) ...
                ./ (w(j) * denom);
            Vref(j + 1, mask) = rho(mask) .* (reduced_flux(mask) - xi(j)) ...
                ./ (w(j + 1) * denom);
            assigned(mask) = true;
        end
    end

    Vref = max(Vref, 0);
    feasible = feasible | assigned;
end

function [xi, w, R, L] = pn_basis(N)
    persistent cache

    if isempty(cache) || cache.N ~= N
        [xi, w] = gausslegendre(N + 1);
        P = legpoly_eval(xi, N);
        G = diag(2 ./ (2 .* (0:N) + 1));
        R = P * diag(w);
        L = (G \ P).';
        cache = struct('N', N, 'xi', xi, 'w', w, 'R', R, 'L', L);
    else
        xi = cache.xi;
        w = cache.w;
        R = cache.R;
        L = cache.L;
    end
end
