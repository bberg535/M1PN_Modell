function [u, uPN] = simulate_periodic_m1pn(rho0, dt, dz, num_steps, NPN, ilim, sigma_a, sigma_s, source_strength)
%SIMULATE_PERIODIC_M1PN Advance the coupled periodic M1PN model.

    nec = numel(rho0);
    nne = nec + 1;

    sigma_a = extend_periodic_field(sigma_a, nec);
    sigma_s = extend_periodic_field(sigma_s, nec);
    flux = @(state) ([state(2,:); state(3,:)]);

    u = zeros(2, nne);
    uPN = zeros(NPN + 1, nne);
    u(1, 1:nec) = reshape(rho0, 1, []);
    uPN(1, 1:nec) = reshape(rho0, 1, []);
    u(:, nne) = u(:, 1);
    uPN(:, nne) = uPN(:, 1);

    for step = 1:num_steps
        step_info = struct('mode', 'M1PN', 'step', step, 'time', step * dt);
        [u, uPN] = heuns_method(u, uPN, dt, dz, flux, nec, ...
            sigma_a, sigma_s, ilim, NPN, source_strength, step_info);
    end
end

function field = extend_periodic_field(field, nec)
    if isscalar(field)
        field = field * ones(1, nec);
    else
        field = reshape(field, 1, []);
    end

    if numel(field) ~= nec
        error('simulate_periodic_m1pn:dimensionMismatch', ...
            'Expected field length %d.', nec);
    end

    field = [field, field(1)];
end
