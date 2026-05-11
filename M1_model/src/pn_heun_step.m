function uPN = pn_heun_step(uPN, dt, dz, N, method, sigma_a, sigma_s, q, boundary, psi_boundary)
%PN_HEUN_STEP Advance the PN auxiliary state by one explicit Heun step.

    if nargin < 9 || isempty(boundary)
        boundary = 'periodic';
    end
    if nargin < 10 || isempty(psi_boundary)
        psi_boundary = 0.0;
    end

    k1 = transport_rhs(uPN, N, method, dt, dz, boundary, psi_boundary);
    uPN1 = uPN + dt * k1;
    uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, q);

    k2 = transport_rhs(uPN1, N, method, dt, dz, boundary, psi_boundary);
    uPN2 = uPN1 + dt * k2;
    uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, q);

    uPN = 0.5 * (uPN + uPN2);
end

function rhs = transport_rhs(uPN, N, method, dt, dz, boundary, psi_boundary)
    flux = pn_interface_flux(uPN, N, method, dt, dz, boundary, psi_boundary);
    rhs = -(flux(:, 2:end) - flux(:, 1:end-1)) / dz;
end
