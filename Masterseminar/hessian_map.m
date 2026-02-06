function H = hessian_map(alpha, Pfull, w)
% H = < b b^T exp(alpha·b) >
z  = exp( (alpha.' * Pfull).');  % (nq x 1)
ws = (w .* z).';                 % (1 x nq)
Bws = Pfull .* ws;               % (nM x nq)
H = Bws * Pfull.';               % (nM x nM)
end