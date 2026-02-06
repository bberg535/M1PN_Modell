function psi_new = euler_advection(psi, mu, dt, dx, psi_vac)
% Upwind FV for each mu
[Nmu, Nx] = size(psi);
psi_new = psi;

% Boundary incoming vacuum
psiL_in = psi_vac;
psiR_in = psi_vac;

for j=1:Nmu
    a = mu(j);
    if a >= 0
        % flux at i+1/2 = a*psi(i)
        % left boundary uses incoming psiL_in
        flux = a * [psiL_in, psi(j,:)];           % length Nx+1 (interfaces)
        psi_new(j,:) = psi(j,:) - (dt/dx) * (flux(2:end) - flux(1:end-1));
    else
        % flux at i+1/2 = a*psi(i+1)
        flux = a * [psi(j,:), psiR_in];           % length Nx+1
        psi_new(j,:) = psi(j,:) - (dt/dx) * (flux(2:end) - flux(1:end-1));
    end
end
end
