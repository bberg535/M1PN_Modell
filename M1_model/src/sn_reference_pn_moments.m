function moments = sn_reference_pn_moments(ref, max_degree)
%SN_REFERENCE_PN_MOMENTS Project an S_N reference onto the PN moment basis.

    if nargin < 2
        max_degree = 2;
    end
    if ~isstruct(ref) || ~isfield(ref, 'psi') || ~isfield(ref, 'mu') || ~isfield(ref, 'w')
        error('sn_reference_pn_moments:invalidReference', ...
            'Expected ref to contain psi, mu and w.');
    end
    if ~isscalar(max_degree) || max_degree < 0 || abs(max_degree - round(max_degree)) > 0
        error('sn_reference_pn_moments:invalidDegree', ...
            'Expected max_degree to be a non-negative integer.');
    end

    P = legpoly_eval(ref.mu, max_degree);
    weighted_psi = reshape(ref.w, [], 1) .* ref.psi;
    moments = P * weighted_psi;
end
