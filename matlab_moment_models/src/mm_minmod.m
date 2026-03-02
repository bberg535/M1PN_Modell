function y = mm_minmod(a, b, c)
%MM_MINMOD Component-wise minmod for vectors of same size.

sameSign = (sign(a) == sign(b)) & (sign(b) == sign(c));
mag = min(abs(a), min(abs(b), abs(c)));
y = zeros(size(a));
y(sameSign) = sign(a(sameSign)) .* mag(sameSign);

end
