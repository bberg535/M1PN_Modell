function [Fpos, Fneg] = flux_parts_MN(alpha, quad, Ppos, Pneg, wmu_pos, wmu_neg)
% Works with alpha either (nM x Nx) or (nM x 1)
if isvector(alpha), alpha = alpha(:); end
Nx = size(alpha,2);

psi_pos = exp( (alpha.' * Ppos).');     % (nq_pos x Nx)
psi_neg = exp( (alpha.' * Pneg).');     % (nq_neg x Nx)

Fpos = Ppos * (wmu_pos .* psi_pos);     % (nM x Nx)
Fneg = Pneg * (wmu_neg .* psi_neg);     % (nM x Nx)
end
