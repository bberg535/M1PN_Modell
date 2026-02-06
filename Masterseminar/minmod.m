function mm = minmod(a,b,c)
mm = zeros(size(a));

mask = (sign(a)==sign(b)) & (sign(a)==sign(c));   % elementwise mask
mm(mask) = sign(a(mask)) .* min(abs(a(mask)), min(abs(b(mask)), abs(c(mask))));

end