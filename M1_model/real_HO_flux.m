function f = real_HO_flux(uPN, u, dt, dx, method, f, fj, Nz)
    %lambda = max(abs(fj(ulc)),abs(fj(urc)));
    lambda = 1;

    HO_Moments = uPN(1:2,:);

    HO_Flux = method(HO_Moments, dt, dx, f, lambda, uPN(3,:), Nz);
    %HO_psi2 = uPN(3,:) ./ max(uPN(1,:),1e-10);

    %HO_Flux = method(u, dt, dx, f, lambda, HO_psi2, Nz);

    %LO_psi2 = calc_psi2(u);

    HO_psi2 = (uPN(1,:) + 2 .* uPN(3,:)) / 3;

    LO_Flux = method(u, dt, dx, f, lambda, HO_psi2, Nz);

    f = calc_g_star(HO_Flux, LO_Flux, f, fj, Nz, u, HO_psi2);
end