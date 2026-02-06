function rho = solve_moment_model(opts)
N   = opts.N;
nM  = N+1;
Nx  = opts.Nx;
dx  = opts.dx;
dt  = opts.dt;
Nt  = opts.Nt;

% Basis evaluation on quadrature nodes for kinetic flux
quad = opts.quad;
Ppos = legendreP_matrix(N, quad.mu_pos);   % (nM x nq_pos)
Pneg = legendreP_matrix(N, quad.mu_neg);   % (nM x nq_neg)
Pfull= legendreP_matrix(N, quad.mu_full);  % (nM x nq_full)

wmu_pos = quad.w_pos(:) .* quad.mu_pos(:);
wmu_neg = quad.w_neg(:) .* quad.mu_neg(:);

% Initial moments: isotropic vacuum everywhere + Dirac distributed to 2 cells (paper)
u = zeros(nM, Nx);
u(1,:) = opts.rho_vac;

midL = Nx/2;
midR = Nx/2 + 1;
u(1, [midL midR]) = opts.rho_vac + 1/dx;  % since rho = 2*(1/(2dx)) = 1/dx

% For MN: store multipliers alpha per cell as cache (initial isotropic guess)
alpha_cache = zeros(nM, Nx);
alpha_cache(1,:) = log(max(u(1,:)/2, 1e-300));  % isotropic alpha0 = log(rho/2)

% Strang Splitting
for it = 1:Nt
    % u_t = s half step (4.1b)
    u = source_halfstep(u, dt/2, opts.sigma_s, opts.sigma_a, opts.Q0);

    % u_t + d_x f_1 + d_y f_2 + d_z f3 = 0 with Heuns method
    [u1, alpha_cache] = euler_transport(u, alpha_cache, dt, dx, opts, Ppos, Pneg, Pfull, wmu_pos, wmu_neg);
    [u2, alpha_cache] = euler_transport(u1, alpha_cache, dt, dx, opts, Ppos, Pneg, Pfull, wmu_pos, wmu_neg);
    u  = 0.5*u + 0.5*u2;

    % u_t = s half step
    u = source_halfstep(u, dt/2, opts.sigma_s, opts.sigma_a, opts.Q0);
end

rho = u(1,:); % for full moments, rho = u0
end