function flux = m1_physical_flux(u)
%M1_PHYSICAL_FLUX Evaluate the M1 physical flux using the Eddington closure.

    flux = [u(2, :); calc_psi2(u)];
end
