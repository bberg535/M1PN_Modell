function flux = m1_interface_flux(u, method, dt, dz, boundary, psi_boundary)
%M1_INTERFACE_FLUX Build interface fluxes for the M1 system.
% For the active coupling this serves as the low-order closed M1 transport
% operator; the HO correction is added only afterwards in real_HO_flux.

    if nargin < 5 || isempty(boundary)
        boundary = 'periodic';
    end
    if nargin < 6 || isempty(psi_boundary)
        psi_boundary = 0.0;
    end

    nec = size(u, 2);
    if nec < 1
        error('m1_interface_flux:emptyState', 'Expected at least one spatial cell.');
    end

    [u_left, u_right, cell_left_idx, cell_right_idx, is_boundary] = ...
        build_interface_states(u, boundary, psi_boundary);
    flux_left = m1_physical_flux(u_left);
    flux_right = m1_physical_flux(u_right);

    fCD = 0.5 * (flux_right + flux_left);
    fLF = fCD - 0.5 * (u_right - u_left);
    fLW = fCD - 0.5 * (dt / dz) * (u_right - u_left);

    if method == -1
        flux = fLF;
        return;
    end

    if method == -2
        flux = fLW;
        return;
    end

    flux = fLF;
    fAe = fLF - fLW;

    [state_min, state_max] = local_state_bounds(u, boundary, psi_boundary);
    for k = 1:size(flux, 2)
        if is_boundary(k)
            continue;
        end

        left_idx = cell_left_idx(k);
        right_idx = cell_right_idx(k);
        wbar = 0.5 * (u_right(:, k) + u_left(:, k)) - 0.5 * ...
            (flux_right(:, k) - flux_left(:, k));
        left_max = state_max(:, left_idx);
        left_min = state_min(:, left_idx);
        right_max = state_max(:, right_idx);
        right_min = state_min(:, right_idx);

        correction = fAe(:, k);
        correction = min(max(0, correction), min(left_max - wbar, wbar - right_min)) ...
            + max(min(0, correction), max(left_min - wbar, wbar - right_max));
        flux(:, k) = fLF(:, k) - correction;
    end

    if strcmpi(boundary, 'periodic')
        flux(:, end) = flux(:, 1);
    end
end

function [u_left, u_right, left_idx, right_idx, is_boundary] = build_interface_states(u, boundary, psi_boundary)
    nec = size(u, 2);
    vacuum_state = [2 * psi_boundary; 0.0];

    switch lower(boundary)
        case 'periodic'
            u_left = [u(:, end), u];
            u_right = [u, u(:, 1)];
            left_idx = [nec, 1:nec];
            right_idx = [1:nec, 1];
            is_boundary = false(1, nec + 1);

        case 'vacuum'
            u_left = [vacuum_state, u];
            u_right = [u, vacuum_state];
            left_idx = [0, 1:nec];
            right_idx = [1:nec, 0];
            is_boundary = false(1, nec + 1);
            is_boundary([1, nec + 1]) = true;

        otherwise
            error('m1_interface_flux:unsupportedBoundary', ...
                'Unsupported boundary condition "%s".', boundary);
    end
end

function [state_min, state_max] = local_state_bounds(u, boundary, psi_boundary)
    nec = size(u, 2);
    vacuum_state = [2 * psi_boundary; 0.0];

    switch lower(boundary)
        case 'periodic'
            prev = [u(:, end), u(:, 1:nec-1)];
            next = [u(:, 2:nec), u(:, 1)];

        case 'vacuum'
            prev = [vacuum_state, u(:, 1:nec-1)];
            next = [u(:, 2:nec), vacuum_state];

        otherwise
            error('m1_interface_flux:unsupportedBoundary', ...
                'Unsupported boundary condition "%s".', boundary);
    end

    state_min = min(u, min(prev, next));
    state_max = max(u, max(prev, next));
end
