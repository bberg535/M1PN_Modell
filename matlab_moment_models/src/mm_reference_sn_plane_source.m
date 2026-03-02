function ref = mm_reference_sn_plane_source(cfg_ref)
%MM_REFERENCE_SN_PLANE_SOURCE High-resolution S_N reference for plane source benchmark.
% Note:
% - Seminarquelle 2 compares against an analytical reference in Section 6.1.1.
% - Here we approximate that reference numerically via a fine S_N solve.

zL = cfg_ref.domain(1);
zR = cfg_ref.domain(2);
tf = cfg_ref.tf;
nCells = cfg_ref.n_cells;
nMu = cfg_ref.n_mu;

sigma_s = cfg_ref.sigma_s;
sigma_a = cfg_ref.sigma_a;
Q = cfg_ref.Q;
psi_vac = cfg_ref.psi_vac_density;

edges = linspace(zL, zR, nCells + 1);
z = 0.5 * (edges(1:end-1) + edges(2:end));
dz = edges(2) - edges(1);

[mu, w] = gauss_legendre(nMu, -1, 1);
maxSpeed = max(abs(mu));

psi = psi_vac * ones(nMu, nCells);

% Plane-source initialization with split Dirac mass (Section 5.5 / 6.1.1).
[~, iCenter] = min(abs(z));
if z(iCenter) <= 0
    iL = iCenter;
    iR = min(nCells, iCenter + 1);
else
    iR = iCenter;
    iL = max(1, iCenter - 1);
end
psi(:, iL) = psi(:, iL) + 1.0 / (2.0 * dz);
psi(:, iR) = psi(:, iR) + 1.0 / (2.0 * dz);

t = 0.0;
while t < tf - 1e-14
    dt = min(cfg_ref.cfl * dz / maxSpeed, tf - t);

    k1 = rhs_sn(psi);
    psi1 = max(psi + dt * k1, 0.0);

    k2 = rhs_sn(psi1);
    psi = max(psi + 0.5 * dt * (k1 + k2), 0.0);

    t = t + dt;
end

rho = w.' * psi;

ref = struct();
ref.z = z(:);
ref.rho = rho(:);
ref.t = t;
ref.mu = mu;
ref.w = w;

    function rhs = rhs_sn(psi_state)
        rho_state = w.' * psi_state;
        % Slab transport source/collision operator from Eq. (2.13) with isotropic scattering.
        src = sigma_s * (0.5 * rho_state - psi_state) - sigma_a * psi_state + Q;

        adv = zeros(size(psi_state));
        for m = 1:nMu
            mmu = mu(m);
            row = psi_state(m, :);

            if mmu >= 0
                F = zeros(1, nCells + 1);
                F(1) = mmu * psi_vac;
                F(2:end) = mmu * row;
            else
                F = zeros(1, nCells + 1);
                F(1:nCells) = mmu * row;
                F(end) = mmu * psi_vac;
            end

            adv(m, :) = -(F(2:end) - F(1:end-1)) / dz;
        end

        rhs = adv + src;
    end

end

function [x, w] = gauss_legendre(n, a, b)
i = (1:n-1)';
beta = i ./ sqrt(4 * i.^2 - 1);
T = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(T);
x = diag(D);
[x, idx] = sort(x);
V = V(:, idx);
w = 2 * (V(1, :)').^2;

x = 0.5 * ((b - a) * x + (a + b));
w = 0.5 * (b - a) * w;
end
