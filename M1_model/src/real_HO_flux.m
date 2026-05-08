function [flux_star, limiter_diag, flux_ho, flux_lo] = real_HO_flux(uPN, u, N, dt, dz, pn_ho_method, boundary, psi_boundary)
%REAL_HO_FLUX Build the limited HOLO flux for the coupled M1PN model.
% Low-order is always the closed M1 LLF flux. The PN system only provides
% the high-order candidate flux before fail-safe limiting.

    if size(u, 2) ~= size(uPN, 2)
        error('real_HO_flux:dimensionMismatch', ...
            'Expected u and uPN to use the same number of spatial cells.');
    end

    flux_lo = m1_interface_flux(u, -1, dt, dz, boundary, psi_boundary);
    pn_flux = pn_interface_flux(uPN, N, pn_ho_method, dt, dz, boundary, psi_boundary);
    flux_ho = pn_flux(1:2, :);
    [flux_star, limiter_diag] = failsafe_limit_density(u, flux_lo, flux_ho, dt, dz, boundary);
end
