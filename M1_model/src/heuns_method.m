function [u, uPN] = heuns_method(u, uPN, dt, dz, flux, nec, sigma_a, sigma_s, ilim, N, source_strength, step_info)
    nne = nec + 1;
    if size(u, 2) ~= nne || size(uPN, 2) ~= nne
        error('heuns_method:dimensionMismatch', ...
            'Expected u and uPN to have %d columns (nec+1).', nne);
    end
    if ilim ~= -1
        error('heuns_method:unsupportedCoupledMethod', ...
            ['Coupled M1PN currently supports only ilim=-1 (LLF). ' ...
            'Requested ilim=%d (%s).'], ilim, method_name(ilim));
    end
    if nargin < 12 || ~isstruct(step_info)
        step_info = struct();
    end
    if ~isfield(step_info, 'mode')
        step_info.mode = 'M1PN';
    end

    u = enforce_periodic(u, nec);
    uPN = enforce_periodic(uPN, nec);
    source = source_term(nne, source_strength);
    pn_source = source(1, :);

    %% Stufe 1
    rhs1 = build_rhs(uPN, u, dt, dz, ilim, flux, nec, source, sigma_a, sigma_s);
    u1 = u + dt * rhs1;
    u1 = enforce_periodic(u1, nec);

    k1PN = advance_pn(uPN, N, ilim, dt, dz);
    uPN1 = uPN + dt * k1PN;
    uPN1 = pn_relax_isotropic(uPN1, N, dt, sigma_a, sigma_s, pn_source);
    uPN1 = enforce_periodic(uPN1, nec);

    %% Stufe 2
    rhs2 = build_rhs(uPN1, u1, dt, dz, ilim, flux, nec, source, sigma_a, sigma_s);
    u2 = u1 + dt * rhs2;
    u2 = enforce_periodic(u2, nec);

    k2PN = advance_pn(uPN1, N, ilim, dt, dz);
    uPN2 = uPN1 + dt * k2PN;
    uPN2 = pn_relax_isotropic(uPN2, N, dt, sigma_a, sigma_s, pn_source);
    uPN2 = enforce_periodic(uPN2, nec);

    %% Zeitschritt
    u = 0.5 * (u + u2);
    u = enforce_periodic(u, nec);

    uPN = 0.5 * (uPN + uPN2);
    uPN = enforce_periodic(uPN, nec);

    assert_valid_m1pn_state(u, uPN, nec, step_info);
end

function rhs = build_rhs(uPN, u, dt, dz, ilim, flux, nec, source, sigma_a, sigma_s)
    im1c = [nec, 1:nec-1];
    Gstar = real_HO_flux(uPN, u, dt, dz, ilim, flux, nec, [], [], sigma_a, sigma_s);

    rhs = zeros(2, nec + 1);
    react = reaction(u(:, 1:nec), sigma_a(1:nec), sigma_s(1:nec));
    rhs(:, 1:nec) = -(Gstar(:, 1:nec) - Gstar(:, im1c)) / dz + source(:, 1:nec) - react;
    rhs(:, nec + 1) = rhs(:, 1);
end

function u = enforce_periodic(u, nec)
    u(:, nec + 1) = u(:, 1);
end

function name = method_name(ilim)
    switch ilim
        case 1
            name = 'MCL';
        case -1
            name = 'LLF';
        otherwise
            name = 'LW';
    end
end
