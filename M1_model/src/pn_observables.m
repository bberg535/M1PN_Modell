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
end
