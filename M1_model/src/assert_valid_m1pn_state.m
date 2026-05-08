function assert_valid_m1pn_state(u, uPN, nec, meta)
%ASSERT_VALID_M1PN_STATE Validate the coupled M1PN end-of-step state.

    if nargin < 4 || ~isstruct(meta)
        meta = struct();
    end

    mode_name = get_field_or(meta, 'mode', 'M1PN');
    step_idx = get_field_or(meta, 'step', []);
    time_value = get_field_or(meta, 'time', []);

    rho = reshape(u(1, 1:nec), 1, []);
    mom = reshape(u(2, 1:nec), 1, []);
    obs = pn_observables(uPN(:, 1:nec));
    closure_obs = pn_closure_observables(uPN(:, 1:nec), u(:, 1:nec));
    pn_m2 = reshape(closure_obs.m2, 1, []);
    pn_min = reshape(obs.min_reconstruction, 1, []);

    reduced_flux = nan(1, nec);
    flux_mask = abs(rho) > 1e-14;
    reduced_flux(flux_mask) = mom(flux_mask) ./ rho(flux_mask);

    [reason, idx] = find_failure(u(:, 1:nec), uPN(:, 1:nec), rho, mom, obs);
    if isempty(reason)
        return;
    end

    context_parts = {};
    if ~isempty(step_idx)
        context_parts{end + 1} = sprintf('step %d', step_idx);
    end
    if ~isempty(time_value)
        context_parts{end + 1} = sprintf('t = %.6g', time_value);
    end

    if isempty(context_parts)
        context_str = '';
    else
        context_str = sprintf(' (%s)', strjoin(context_parts, ', '));
    end

    error('assert_valid_m1pn_state:invalidState', ...
        ['%s invalid final state%s, cell %d: reason=%s, rho=%.16e, ' ...
        'j=%.16e, j/rho=%.16e, pn_m2=%.16e, min PN reconstruction=%.16e'], ...
        mode_name, context_str, idx, reason, rho(idx), mom(idx), ...
        reduced_flux(idx), pn_m2(idx), pn_min(idx));
end

function [reason, idx] = find_failure(u, uPN, rho, mom, obs)
    tol = 1e-12;
    reason = '';
    idx = [];

    invalid_u = ~isfinite(u) | ~isreal(u);
    if any(invalid_u(:))
        [~, idx] = first_spatial_index(invalid_u);
        reason = 'M1 contains non-finite or complex values';
        return;
    end

    invalid_uPN = ~isfinite(uPN) | ~isreal(uPN);
    if any(invalid_uPN(:))
        [~, idx] = first_spatial_index(invalid_uPN);
        reason = 'PN contains non-finite or complex values';
        return;
    end

    neg_rho = rho < -tol;
    if any(neg_rho)
        [~, idx] = min(rho);
        reason = 'negative density';
        return;
    end

    realizability_excess = abs(mom) - rho;
    if any(realizability_excess > tol)
        [~, idx] = max(realizability_excess);
        reason = '|j| exceeds rho';
        return;
    end

    invalid_m2 = ~isfinite(obs.m2) | ~isreal(obs.m2);
    if any(invalid_m2)
        idx = find(invalid_m2, 1, 'first');
        reason = 'PN-based m2 is non-finite or complex';
        return;
    end

    invalid_recon = ~isfinite(obs.angular_values) | ~isreal(obs.angular_values);
    if any(invalid_recon(:))
        [~, idx] = first_spatial_index(invalid_recon);
        reason = 'PN reconstruction is non-finite or complex';
    end
end

function [row_idx, col_idx] = first_spatial_index(mask)
    [row_idx, col_idx] = find(mask, 1, 'first');
end

function value = get_field_or(s, field_name, default_value)
    if isfield(s, field_name)
        value = s.(field_name);
    else
        value = default_value;
    end
end
