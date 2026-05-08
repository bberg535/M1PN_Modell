function flux_u = pn_interface_flux(uPN, N, method, dt, dz, boundary, psi_boundary)
%PN_INTERFACE_FLUX Build interface fluxes for the PN transport system.
% The transport step is assembled in the diagonalized angular variables v.
% For the HOLO coupling this routine is the single source of PN interface
% fluxes, both for pure PN runs and for PN-based HO candidate fluxes.

    if nargin < 6 || isempty(boundary)
        boundary = 'periodic';
    end
    if nargin < 7 || isempty(psi_boundary)
        psi_boundary = 0.0;
    end

    nec = size(uPN, 2);
    [xi, ~, R, L] = pn_basis(N);
    v = L * uPN;

    [v_left, v_right, left_idx, right_idx, is_boundary, v_left_bc, v_right_bc] = ...
        build_interface_states(v, xi, boundary, psi_boundary);
    flux_left = xi .* v_left;
    flux_right = xi .* v_right;
    lambda = abs(xi);

    fCD = 0.5 * (flux_right + flux_left);
    fLF = fCD - 0.5 * (lambda .* (v_right - v_left));
    fLW = fCD - 0.5 * (dt / dz) * ((lambda .^ 2) .* (v_right - v_left));

    if method == -1
        flux_v = fLF;
    elseif method == -2
        flux_v = fLW;
    else
        flux_v = fLF;
        fAe = fLF - fLW;
        [v_min, v_max] = local_state_bounds(v, boundary, v_left_bc, v_right_bc);
        for k = 1:size(flux_v, 2)
            if is_boundary(k)
                continue;
            end

            left_cell = left_idx(k);
            right_cell = right_idx(k);
            wbar = 0.5 * (v_right(:, k) + v_left(:, k)) .* lambda ...
                - 0.5 * (flux_right(:, k) - flux_left(:, k));
            left_max = v_max(:, left_cell);
            left_min = v_min(:, left_cell);
            right_max = v_max(:, right_cell);
            right_min = v_min(:, right_cell);

            correction = fAe(:, k);
            correction = min(max(0, correction), ...
                min(lambda .* left_max - wbar, wbar - lambda .* right_min)) ...
                + max(min(0, correction), ...
                max(lambda .* left_min - wbar, wbar - lambda .* right_max));
            flux_v(:, k) = fLF(:, k) - correction;
        end
    end

    if strcmpi(boundary, 'periodic')
        flux_v(:, end) = flux_v(:, 1);
    end

    flux_u = R * flux_v;
end

function [v_left, v_right, left_idx, right_idx, is_boundary, ghost_left, ghost_right] = ...
    build_interface_states(v, xi, boundary, psi_boundary)
    nec = size(v, 2);
    xi = xi(:);

    switch lower(boundary)
        case 'periodic'
            ghost_left = v(:, end);
            ghost_right = v(:, 1);
            v_left = [v(:, end), v];
            v_right = [v, v(:, 1)];
            left_idx = [nec, 1:nec];
            right_idx = [1:nec, 1];
            is_boundary = false(1, nec + 1);

        case 'vacuum'
            ghost_left = v(:, 1);
            ghost_left(xi > 0) = psi_boundary;

            ghost_right = v(:, end);
            ghost_right(xi < 0) = psi_boundary;

            v_left = [ghost_left, v];
            v_right = [v, ghost_right];
            left_idx = [0, 1:nec];
            right_idx = [1:nec, 0];
            is_boundary = false(1, nec + 1);
            is_boundary([1, nec + 1]) = true;

        otherwise
            error('pn_interface_flux:unsupportedBoundary', ...
                'Unsupported boundary condition "%s".', boundary);
    end
end

function [v_min, v_max] = local_state_bounds(v, boundary, ghost_left, ghost_right)
    nec = size(v, 2);

    switch lower(boundary)
        case 'periodic'
            prev = [v(:, end), v(:, 1:nec-1)];
            next = [v(:, 2:nec), v(:, 1)];

        case 'vacuum'
            prev = [ghost_left, v(:, 1:nec-1)];
            next = [v(:, 2:nec), ghost_right];

        otherwise
            error('pn_interface_flux:unsupportedBoundary', ...
                'Unsupported boundary condition "%s".', boundary);
    end

    v_min = min(v, min(prev, next));
    v_max = max(v, max(prev, next));
end
