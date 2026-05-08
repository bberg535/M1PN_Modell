function [xi, w, R, L] = pn_basis(N)
%PN_BASIS Return cached quadrature and transformation operators for PN order N.
% The cache keeps all active PN helper routines on the same basis data.

    persistent cache

    if isempty(cache)
        cache = struct('N', {}, 'xi', {}, 'w', {}, 'R', {}, 'L', {});
    end

    match_idx = find([cache.N] == N, 1, 'first');
    if isempty(match_idx)
        [xi, w] = gausslegendre(N + 1);
        P = legpoly_eval(xi, N);
        G = diag(2 ./ (2 .* (0:N) + 1));
        R = P * diag(w);
        L = (G \ P).';
        cache(end + 1) = struct('N', N, 'xi', xi, 'w', w, 'R', R, 'L', L); %#ok<AGROW>
    else
        xi = cache(match_idx).xi;
        w = cache(match_idx).w;
        R = cache(match_idx).R;
        L = cache(match_idx).L;
    end
end
