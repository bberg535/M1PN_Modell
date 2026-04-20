function m1_fv_mcl
    % Wahl des Verfahrens:
    %  ilim =  1  -> MCL
    %  ilim = -1  -> LLF
    %  sonst      -> LW
    ilim   = 1;

    % Diskretisierung
    dt     = 0.0001;
    dz     = 0.0002;        % CFL = dt/dz = 0.5
    Nz     = 2/dz;
    Tfinal = 1.0;

    z      = (dz/2:dz:2-dz/2);

    % Zustand U = [psi0; psi1] (2 x Nz)
    U      = zeros(2, Nz);
    rho0   = [1/dz; zeros(Nz-1,1)];
    U(1,:) = rho0.';        % psi0
    U(2,:) = 0;             % psi1

    % Flux & Wellenzahl-Funktion
    flux  = @m1_flux;       % siehe Funktion unten
    fluxj = @(U) 1;         % max. Wellengeschw. ~ 1

    % Auswahl der FV-Methode
    switch ilim
        case 1
            method = @MCL_rhs;
        case -1
            method = @LLF_rhs;
        otherwise
            method = @LW_rhs;
    end

    % Anfangsmasse
    mass0 = sum(U(1,:)) * dz;

    % Plot
    figure; 
    h = plot(z, U(1,:), z, U(2,:), z, compute_psi2(U));
    legend('\psi_0','\psi_1','\psi_2');
    xlabel('z'); grid on;

    t    = 0.0;
    step = 0;

    while t < Tfinal - 1e-12
        if t + dt > Tfinal
            dt = Tfinal - t;
        end

        step = step + 1;

        % Heun-Schritt auf (psi0, psi1)
        U = heun_step(U, dt, dz, flux, fluxj, Nz, method);

        % Projektion in die Realisierbarkeitsmenge:
        %  psi0 > 0, |psi1| <= psi0
        U = project_realizable(U);

        t = t + dt;

        if mod(step,50) == 0
            psi2 = compute_psi2(U);
            mass = sum(U(1,:)) * dz;
            fprintf('t = %.3f, mass = %.10e, min(psi0) = %.3e\n', ...
                    t, mass, min(U(1,:)));

            set(h(1), 'YData', U(1,:));
            set(h(2), 'YData', U(2,:));
            set(h(3), 'YData', psi2);
            drawnow;
        end
    end
end


%======================================================================
% Heun (SSP-RK2) für U_t = RHS(U)
%======================================================================
function Unew = heun_step(U, dt, dx, flux, fluxj, Nz, method)
    % max. Wellenzahl (hier konstant 1)
    lambda = 1; %#ok<NASGU> % bleibt drin, falls du es später variabel machst

    k1   = method(U, dt, dx, flux, Nz);
    U1   = U + dt * k1;
    k2   = method(U1, dt, dx, flux, Nz);
    Unew = U + 0.5 * dt * (k1 + k2);
end


%======================================================================
% Physikalischer M1-Flux F(U) = [psi1; psi2(psi0,psi1)]
%======================================================================
function F = m1_flux(U)
    psi0 = U(1,:);
    psi1 = U(2,:);

    eps_rho   = 1e-14;
    psi0_safe = max(psi0, eps_rho);

    f   = abs(psi1) ./ psi0_safe;
    % Hier kein hartes f-Clipping: Realisierbarkeit wird über
    % project_realizable(U) sichergestellt.
    chi = (3 + 4*f.^2) ./ (5 + 2*sqrt(4 - 3*f.^2));
    psi2 = chi .* psi0_safe;

    F = zeros(size(U));
    F(1,:) = psi1;
    F(2,:) = psi2;
end


%======================================================================
% Hilfsfunktion nur für den Plot: psi2 aus aktuellem U
%======================================================================
function psi2 = compute_psi2(U)
    psi0 = U(1,:);
    psi1 = U(2,:);

    eps_rho   = 1e-14;
    psi0_safe = max(psi0, eps_rho);
    f   = abs(psi1) ./ psi0_safe;
    chi = (3 + 4*f.^2) ./ (5 + 2*sqrt(4 - 3*f.^2));
    psi2 = chi .* psi0_safe;
end


%======================================================================
% LLF-Rechte Seite: dudt = -(F_{j+1/2} - F_{j-1/2}) / dx
%======================================================================
function dudt = LLF_rhs(U, dt, dx, flux, Nz)
    ip0c = 1:Nz;
    ip1c = [2:Nz 1];
    im1c = [Nz 1:Nz-1];

    ulc = U(:,ip0c);
    urc = U(:,ip1c);

    lambda = 1;

    fC  = 0.5 * (flux(ulc) + flux(urc));
    fLF = fC - 0.5 * lambda * (urc - ulc);

    Fp   = fLF;
    Fm   = fLF(:,im1c);
    dudt = -(Fp - Fm) / dx;
end


%======================================================================
% Lax–Wendroff RHS (einfach, nicht HOLO-korrekt, aber ok zum Testen)
%======================================================================
function dudt = LW_rhs(U, dt, dx, flux, Nz)
    ip0c = 1:Nz;
    ip1c = [2:Nz 1];
    im1c = [Nz 1:Nz-1];

    ulc = U(:,ip0c);
    urc = U(:,ip1c);

    lambda = 1;
    cfl    = dt/dx;

    fC  = 0.5 * (flux(ulc) + flux(urc));
    fLW = fC - 0.5 * cfl * lambda^2 * (urc - ulc);

    Fp   = fLW;
    Fm   = fLW(:,im1c);
    dudt = -(Fp - Fm) / dx;
end


%======================================================================
% MCL-RHS: LLF + begrenzter antidiffusiver Flux
%======================================================================
function dudt = MCL_rhs(U, dt, dx, flux, Nz)
    ip0c = 1:Nz;
    ip1c = [2:Nz 1];
    im1c = [Nz 1:Nz-1];

    ulc = U(:,ip0c);
    urc = U(:,ip1c);

    lambda = 1;
    cfl    = dt/dx;

    lambdaMCL = lambda * ones(size(ulc));

    fC  = 0.5 * (flux(ulc) + flux(urc));
    fLF = fC - 0.5 * lambda * (urc - ulc);
    fLW = fC - 0.5 * cfl * lambda^2 * (urc - ulc);

    fAe = fLF - fLW;

    umax = max(U(:,im1c), max(U(:,ip0c), U(:,ip1c)));
    umin = min(U(:,im1c), min(U(:,ip0c), U(:,ip1c)));

    wbar = 0.5*(urc+ulc).*lambdaMCL - 0.5*(flux(urc)-flux(ulc));

    fAe_pos = min( max(0, fAe), ...
                   min(lambdaMCL.*umax(:,ip0c)-wbar, ...
                       wbar - lambdaMCL.*umin(:,ip1c)) );

    fAe_neg = max( min(0, fAe), ...
                   max(lambdaMCL.*umin(:,ip0c)-wbar, ...
                       wbar - lambdaMCL.*umax(:,ip1c)) );

    fMCL = fLF - (fAe_pos + fAe_neg);

    Fp   = fMCL;
    Fm   = fMCL(:,im1c);
    dudt = -(Fp - Fm) / dx;
end


%======================================================================
% Projektion auf Realisierbarkeitsmenge:
%   psi0 > 0, |psi1| <= psi0
%======================================================================
function Uproj = project_realizable(U)
    psi0 = U(1,:);
    psi1 = U(2,:);

    eps_rho = 1e-14;

    % psi0 minimal positiv
    psi0(psi0 < eps_rho) = eps_rho;

    % |psi1| <= (1 - eps) * psi0
    mask = abs(psi1) > (1 - 1e-10).*psi0;
    psi1(mask) = (1 - 1e-10).*psi0(mask).*sign(psi1(mask));

    Uproj      = U;
    Uproj(1,:) = psi0;
    Uproj(2,:) = psi1;
end
