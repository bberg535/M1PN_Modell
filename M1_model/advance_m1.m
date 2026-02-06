function u = advance_m1(un, dt, dx, f, fj, Nz, sigma_a, sigma_s, method, source_strength)
    % Linker und Rechter Zustand
    ip0c = (1:Nz)'; ip1c = [2:Nz 1]'; im1c = [Nz 1:Nz-1];
    ulc = un(:,ip0c); urc = un(:,ip1c);

    % Wellengeschwindigkeit
    lambda = max(abs(fj(ulc)),abs(fj(urc)));

    % Stufe 1
    psi2_0 = calc_psi2(un);
    f1  = method(un, dt, dx, f, lambda, psi2_0, Nz);
    u1  = un + dt .* (-(f1 - f1(:,im1c)) / dx + source_term(Nz, source_strength) - reaction(un, sigma_a, sigma_s));

    % Stufe 2
    psi2_1 = calc_psi2(u1);
    f2  = method(u1, dt, dx, f, lambda, psi2_1, Nz);
    u2 = u1 + dt .* (-(f2 - f2(:,im1c)) / dx + source_term(Nz, source_strength) - reaction(u1, sigma_a, sigma_s));

    % Heun-Update
    u = (un + u2)/2;
end