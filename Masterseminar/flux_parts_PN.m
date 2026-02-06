function [Fpos, Fneg] = flux_parts_PN(u, N, quad, Ppos, Pneg, wmu_pos, wmu_neg)
nM  = N+1;
if isvector(u), u = u(:); end
Nx  = size(u,2);

% PN-Fluss: Wie im PPN-Paper konstruiere A Matrix, die über die Integrale
% der Legendre-Polynome definiert ist
k = (0:N)';
a_scale = (2*k+1)/2;                 % (nM x 1)
A = a_scale .* u;                    % (nM x Nx)

% Evaluate psi at quadrature nodes
psi_pos = (A.' * Ppos).';            % (nq_pos x Nx) since Ppos is (nM x nq_pos)
psi_neg = (A.' * Pneg).';            % (nq_neg x Nx)

% Flux parts
Fpos = Ppos * (wmu_pos .* psi_pos);  % (nM x Nx)
Fneg = Pneg * (wmu_neg .* psi_neg);  % (nM x Nx) (note: mu_neg negative already)
end

