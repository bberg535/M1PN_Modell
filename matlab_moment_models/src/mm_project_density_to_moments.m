function u = mm_project_density_to_moments(psi_handle, model, quad)
%MM_PROJECT_DENSITY_TO_MOMENTS Project density psi(mu) to model moments.
% References:
% - Moment definition u=<b psi>: Seminarquelle 1 Eq. (2.3), Seminarquelle 2 Eq. (2.13) in slab notation.

psi = psi_handle(quad.mu);
psi = psi(:);

if numel(psi) ~= numel(quad.mu)
    error('psi_handle must return one value per quadrature node.');
end

u = quad.B.' * (quad.w .* psi);

end
