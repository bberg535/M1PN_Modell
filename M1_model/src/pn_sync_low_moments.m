function [u_sync, sync_diag] = pn_sync_low_moments(uPN, target_moments)
%PN_SYNC_LOW_MOMENTS Least-squares projection onto prescribed rho/j moments.

    if size(target_moments, 1) ~= 2
        error('pn_sync_low_moments:invalidTargetShape', ...
            'Expected target_moments to have size 2 x Ncells.');
    end
    if size(uPN, 2) ~= size(target_moments, 2)
        error('pn_sync_low_moments:dimensionMismatch', ...
            'Expected PN state and target moments to share the same number of cells.');
    end
    if size(uPN, 1) < 2
        error('pn_sync_low_moments:insufficientOrder', ...
            'At least a P1 auxiliary basis is required to match rho and j.');
    end

    A = low_moment_map(size(uPN, 1) - 1);
    gram = A * A.';

    pre_mismatch = A * uPN - target_moments;
    lagrange = gram \ pre_mismatch;
    correction = A.' * lagrange;
    u_sync = uPN - correction;
    post_mismatch = A * u_sync - target_moments;

    correction_norm = sqrt(sum(correction .^ 2, 1));
    sync_diag = struct();
    sync_diag.max_pre_mismatch = max(abs(pre_mismatch(:)));
    sync_diag.max_post_mismatch = max(abs(post_mismatch(:)));
    sync_diag.max_correction = max(correction_norm);
    sync_diag.mean_correction = mean(correction_norm);
end

function A = low_moment_map(N)
    persistent cache

    if isempty(cache)
        cache = struct('N', {}, 'A', {});
    end

    match_idx = find([cache.N] == N, 1, 'first');
    if isempty(match_idx)
        [xi, w, ~, L] = pn_basis(N);
        A = [w.' * L; (w(:) .* xi(:)).' * L];
        cache(end + 1) = struct('N', N, 'A', A); %#ok<AGROW>
    else
        A = cache(match_idx).A;
    end
end
