function rho = solve_SN_reference(ref, tf)
% Discrete ordinates reference for slab transport with isotropic scattering.
% Uses exact scattering halfsteps and SSP-RK2 advection.

Nx = ref.Nx;
dx = ref.dx;
dt = ref.dt;
Nt = ceil(tf/dt);
dt = tf/Nt;

% Gauss-Legendre quadrature in mu
[mu, w] = gauss_legendre(ref.Nmu, -1, 1);
mu = mu(:); w = w(:);

% Initial psi: vacuum + Dirac on two center cells
psi = ref.psi_vac * ones(ref.Nmu, Nx);
midL = Nx/2; midR = Nx/2 + 1;
psi(:, [midL midR]) = ref.psi_vac + 1/(2*dx);

for it=1:Nt
    psi = scatter_halfstep(psi, w, dt/2, ref.sigma_s, ref.sigma_a, ref.Q0);

    psi1 = euler_advection(psi, mu, dt, dx, ref.psi_vac);
    psi2 = euler_advection(psi1, mu, dt, dx, ref.psi_vac);
    psi  = 0.5*psi + 0.5*psi2;

    psi = scatter_halfstep(psi, w, dt/2, ref.sigma_s, ref.sigma_a, ref.Q0);
end

rho = (w.' * psi); % rho(x) = <psi> = sum_j w_j psi_j(x)
end