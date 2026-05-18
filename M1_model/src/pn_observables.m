function obs = pn_observables(uPN)
%PN_OBSERVABLES Compute rho, j and m2 from PN moments via quadrature reconstruction.

    N = size(uPN, 1) - 1;
    [xi, w, ~, L] = pn_basis(N);
    v = L * uPN;

    obs = struct();
    obs.rho = (w.' * v).';
    obs.j = ((w(:) .* xi(:)).' * v).';
    obs.m2 = ((w(:) .* (xi(:) .^ 2)).' * v).';
    obs.angular_values = v;
    obs.min_reconstruction = min(v, [], 1).';

    tol = 1.0e-11 * max(1.0, max(abs(uPN(:))));
    rho_err = max(abs(obs.rho - uPN(1, :).'));
    if rho_err > tol
        error('pn_observables:rhoMismatch', ...
            ['u_0 must equal the angular density integral. ' ...
             'Observed max |rho - u_0| = %.6e.'], rho_err);
    end

    if N >= 1
        j_err = max(abs(obs.j - uPN(2, :).'));
        if j_err > tol
            error('pn_observables:jMismatch', ...
                ['u_1 must equal the angular flux integral. ' ...
                 'Observed max |j - u_1| = %.6e.'], j_err);
        end
    end
end
