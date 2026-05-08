function P = evaluate_legendre_polynomial(k, mu)
% Wertet Legendre-Polynom P_k(μ) aus mit Drei-Term-Rekursion
%
% P_0(μ) = 1
% P_1(μ) = μ
% (n+1)*P_{n+1}(μ) = (2n+1)*μ*P_n(μ) - n*P_{n-1}(μ)

if k == 0
    P = ones(size(mu));
    return;
elseif k == 1
    P = mu;
    return;
end

% Rekursion für k ≥ 2
P_km2 = ones(size(mu));   % P_0
P_km1 = mu;               % P_1

for n = 1:(k-1)
    P = ((2*n + 1) * mu .* P_km1 - n * P_km2) / (n + 1);
    P_km2 = P_km1;
    P_km1 = P;
end
end