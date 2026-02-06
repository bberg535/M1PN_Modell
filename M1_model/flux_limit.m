function rhs = flux_limit(HO_Flux, LO_Flux, Nz, u, psi2, dt, dz, source_strength, sigma_a, sigma_s)
    ip0c=[1:Nz]'; ip1c=[2:Nz 1]'; im1c=[Nz 1:Nz-1]';

    rhs = zeros(2, Nz+1);
    flux = [u(2,:); psi2];
    react = reaction(u, sigma_a, sigma_s);
    source = source_term(Nz+1, source_strength);
    fAe = HO_Flux - LO_Flux; % Gleichung (11)

    % Local bounds
    umax = max(u(:,im1c), max(u(:,ip0c), u(:,ip1c))); % Gleichung (24a) mit Gleichung (22)
    umin = min(u(:,im1c), min(u(:,ip0c), u(:,ip1c))); % Gleichung (24a) mit Gleichung (22) 

    % Scaled bar states
    % Im Buch (3.72) aber mit lambda = 2*d_ij multipliziert
    %wbar = 0.5 * (u(:,ip0c)+u(:,ip1c)) - 0.5 * (flux(:,ip1c) - flux(:,ip0c)) + source(:, ip0c) - react(:, ip0c); 
    wbar = 0.5 * (u(:,ip0c)+u(:,ip1c)) - 0.5 * (flux(:,ip1c) - flux(:,ip0c)); 

    % Flux limiting
    % f_ij max (Unter Gleichung (3.77))
    fAe = min(max(0,fAe), min(umax(:,ip0c) - wbar, wbar - umin(:,ip1c))) ...
        + max(min(0,fAe), max(umin(:,ip0c) - wbar, wbar - umax(:,ip1c))); % Gleichung (27)

    % Flux correction
    Gstar = LO_Flux - fAe;
    rhs(:,ip0c) = 1/dz * (wbar(:,im1c) + Gstar(:,im1c) - u(:,ip0c) + wbar(:,ip1c) - Gstar(:,ip1c) - u(:,ip0c));
end
