function [u, uPN] = heuns_method(u, uPN, dt, dz, flux, fluxj, Nz, sigma_a, sigma_s, ilim, N)
    im1c = [Nz, 1:Nz-1]; cfl = dt / dz;

    switch ilim
        case 1
            method = @MCL;
        case -1
            method = @LLF;
        otherwise
            method = @LW;
    end  

    % ---------- STUFE 1 ----------
    % HO+LO Flux für u^n mit uPN^n
    f1 = real_HO_flux(uPN, u, dt, dz, method, flux, fluxj, Nz);

    % Reaktion & Quelle (M1)
    rhs1 = -(f1 - f1(:,im1c))/dz + source_term(Nz) - reaction(u, sigma_a, sigma_s);

    u1 = u + dt * rhs1;

    % PN: erste Stufe (hier brauchst du eine semi-diskrete pn_rhs)
    k1PN = advance_pn(uPN, N, ilim, dt, dz);   % baust du aus deinem PPN-Operator
    uPN1 = uPN + dt * k1PN;

    % ---------- STUFE 2 ----------
    f2 = real_HO_flux(uPN1, u1, dt, dz, method, flux, fluxj, Nz);
    rhs2 = -(f2 - f2(:,im1c))/dz + source_term(Nz) - reaction(u1, sigma_a, sigma_s);

    u2 = u1 + dt * rhs2;

    % PN zweite Stufe
    k2PN = advance_pn(uPN1, N, ilim, dt, dz);
    uPN2 = uPN1 + dt * k2PN;

    % ---------- Heun-Mittelung ----------
    u   = 0.5 * (u + u2);
    uPN = 0.5 * (uPN + uPN2);
end