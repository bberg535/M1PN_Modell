function [flux_star, limiter] = failsafe_limit_density(u, flux_lo, flux_ho, dt, dz, boundary)
%FAILSAFE_LIMIT_DENSITY Fail-safe flux limiting for the density component.
% This is the V1 realization of the monolithic fail-safe limiter from
% Section 2.5.6.2 of "Property-Preserving Numerical Schemes for
% Conservation Laws", restricted to the invariant-domain condition rho >= 0.

    if nargin < 6 || isempty(boundary)
        boundary = 'periodic';
    end

    nec = size(u, 2);
    cfl = dt / dz;
    % Antidiffusive correction that would upgrade LO to HO without limiting.
    anti_flux = flux_ho - flux_lo;
    center_flux = m1_physical_flux(u);
    rho = reshape(u(1, :), 1, []);
    alpha = ones(1, size(flux_lo, 2));

    for k = 1:numel(alpha)
        [left_idx, right_idx] = adjacent_cells(k, nec, boundary);
        a_rho = anti_flux(1, k);
        alpha_k = 1.0;

        % Limit whichever cell would lose density under the correction.
        if left_idx > 0 && a_rho < 0
            rho_left = rho(left_idx) - 2 * cfl * ...
                (flux_lo(1, k) - center_flux(1, left_idx));
            alpha_k = min(alpha_k, admissible_alpha(rho_left, -a_rho, cfl));
        end

        if right_idx > 0 && a_rho > 0
            rho_right = rho(right_idx) - 2 * cfl * ...
                (center_flux(1, right_idx) - flux_lo(1, k));
            alpha_k = min(alpha_k, admissible_alpha(rho_right, a_rho, cfl));
        end

        alpha(k) = max(0.0, min(1.0, alpha_k));
    end

    if strcmpi(boundary, 'periodic')
        alpha(end) = alpha(1);
    end

    flux_star = flux_lo + anti_flux .* alpha;
    limiter = summarize_limiter(alpha, boundary);
end

function alpha = admissible_alpha(rho_state, a_rho, cfl)
    tol = 1e-14;
    if a_rho <= tol
        alpha = 1.0;
        return;
    end
    alpha = rho_state / (2 * cfl * a_rho);
end

function [left_idx, right_idx] = adjacent_cells(interface_idx, nec, boundary)
    switch lower(boundary)
        case 'periodic'
            left_idx = mod(interface_idx - 2, nec) + 1;
            right_idx = mod(interface_idx - 1, nec) + 1;

        case 'vacuum'
            if interface_idx == 1
                left_idx = 0;
                right_idx = 1;
            elseif interface_idx == nec + 1
                left_idx = nec;
                right_idx = 0;
            else
                left_idx = interface_idx - 1;
                right_idx = interface_idx;
            end

        otherwise
            error('failsafe_limit_density:unsupportedBoundary', ...
                'Unsupported boundary condition "%s".', boundary);
    end
end

function limiter = summarize_limiter(alpha, boundary)
    tol = 1e-12;
    if strcmpi(boundary, 'periodic')
        alpha_unique = alpha(1:end-1);
    else
        alpha_unique = alpha;
    end

    limiter = struct();
    limiter.alpha = alpha;
    limiter.min_alpha = min(alpha_unique);
    limiter.max_alpha = max(alpha_unique);
    limiter.num_limited = nnz(alpha_unique < 1 - tol);
end
