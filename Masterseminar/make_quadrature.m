function quad = make_quadrature(nq)
% Gauss-Legendre on [-1,1] and on [-1,0],[0,1] for kinetic splitting
[mu_full, w_full] = gauss_legendre(nq, -1, 1);
[mu_pos,  w_pos ] = gauss_legendre(nq,  0, 1);
[mu_neg,  w_neg ] = gauss_legendre(nq, -1, 0);

quad.mu_full = mu_full(:); quad.w_full = w_full(:);
quad.mu_pos  = mu_pos(:);  quad.w_pos  = w_pos(:);
quad.mu_neg  = mu_neg(:);  quad.w_neg  = w_neg(:);
end
