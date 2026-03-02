function [flux, flux_state] = mm_kinetic_flux(uL, uR, model, quad, flux_cfg)
%MM_KINETIC_FLUX Kinetic interface flux g(uL,uR) for slab models.
% Reference: Seminarquelle 2 Eq. (4.15), using half-range integrals <·>+ and <·>-.

if nargin < 5
    flux_cfg = struct();
end

opt_cfg = get_field_or(flux_cfg, 'opt_cfg', struct());
cacheL = get_field_or(flux_cfg, 'cache_left', struct());
cacheR = get_field_or(flux_cfg, 'cache_right', struct());

uL = uL(:);
uR = uR(:);

alphaL = [];
alphaR = [];
infoL = struct();
infoR = struct();

if model.needs_entropy
    [alphaL, infoL] = mm_entropy_dual_solve(uL, model, quad, opt_cfg, cacheL);
    [alphaR, infoR] = mm_entropy_dual_solve(uR, model, quad, opt_cfg, cacheR);

    psiL_plus = exp(min(quad.B_plus * alphaL, 700));
    psiR_minus = exp(min(quad.B_minus * alphaR, 700));
else
    alphaL = model.mass_matrix \ uL;
    alphaR = model.mass_matrix \ uR;

    psiL_plus = quad.B_plus * alphaL;
    psiR_minus = quad.B_minus * alphaR;
end

flux = quad.B_plus.' * (quad.w_plus .* quad.mu_plus .* psiL_plus) + ...
       quad.B_minus.' * (quad.w_minus .* quad.mu_minus .* psiR_minus);

flux_state = struct('alpha_left', alphaL, 'alpha_right', alphaR, ...
    'info_left', infoL, 'info_right', infoR);

end

function v = get_field_or(s, name, default)
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
