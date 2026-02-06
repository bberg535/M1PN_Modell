function u = moment_map(alpha, Pfull, w)
% u = < b exp(alpha·b) >
z = exp( (alpha.' * Pfull).');   % (nq x 1)
u = Pfull * (w .* z);            % (nM x 1)
end
